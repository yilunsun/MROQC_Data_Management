library(Hmisc)
library(tidyverse)
library(lubridate) 

setwd(getwd())
setwd('~/Dropbox/GSRA/MROQC')
### read in data
w_rp2 = read_csv("w_rp2.csv");

plotdata = w_rp2 %>% 
  select("RT", "rp2", "w_rp2", "rt_st_Date") %>%
  mutate(Year = year(as.Date(rt_st_Date, "%d%B%Y"))) %>%
  filter(complete.cases(.)) %>%
  mutate(weight = w_rp2*n()/sum(w_rp2)) %>%
  select(-rt_st_Date, - w_rp2)

plotdata_unwt = plotdata %>%
  group_by(Year, RT) %>%
  summarise(Prob = mean(rp2), upper = binconf(sum(rp2),n(),method="wilson")[3], lower = binconf(sum(rp2),n(),method="wilson")[2])

plotdata_wt = plotdata %>%
  group_by(Year, RT) %>%
  summarise(Prob = sum(rp2*weight)/sum(weight), upper = binconf(sum(rp2*weight),sum(weight),method="wilson")[3], lower = binconf(sum(rp2*weight),sum(weight),method="wilson")[2])

### start plot 1
ggplot(data = plotdata_unwt, aes(x=Year, y=Prob, color=RT, group = RT)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=.2,
                position=position_dodge(0.05)) +
  ylab("Probability of Gr 2 or Higher Pneumonitis") +
  theme(axis.title.x = element_text(face="bold", colour="black",size=18),axis.text.x=element_text(size=20,face="bold",colour="black"))+
  theme(axis.title.y = element_text(face="bold", colour="black",size=18),axis.text.y = element_text(face="bold", colour="black",size=20))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size = 16, face = "bold"))+
  theme(legend.key.size=unit(1.2,"cm"))

### start plot 2
ggplot(data = plotdata_wt, aes(x=Year, y=Prob, color=RT, group = RT)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper),width=.2,
                position=position_dodge(0.05)) +
  ylab("Probability of Gr 2 or Higher Pneumonitis") +
  theme(axis.title.x = element_text(face="bold", colour="black",size=18),axis.text.x=element_text(size=20,face="bold",colour="black"))+
  theme(axis.title.y = element_text(face="bold", colour="black",size=18),axis.text.y = element_text(face="bold", colour="black",size=20))+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size = 16, face = "bold"))+
  theme(legend.key.size=unit(1.2,"cm"))