library(tidyverse)

setwd("~/Dropbox/GSRA/MROQC")

weightloss = read_csv("weight_loss.csv")

weightloss = weightloss %>%
  mutate(loss = ifelse(is.na(Weight_L4), -(weight_correct_l8 - Weight_L1)/Weight_L1, -(weight_correct_l8 - Weight_L4)/Weight_L4)) %>%
  filter(!is.na(loss))

psych::describe(weightloss$loss)

ordered = weightloss %>%
  group_by(hosp) %>%
  summarise(average = median(loss)) %>%
  arrange(average) %>%
  select(hosp)

ordered = unlist(ordered)

weightloss$hosp = factor(weightloss$hosp, levels=ordered)

ggplot(data=weightloss) +
  geom_boxplot(aes(x = hosp, y = loss*100)) +
  ylab('% Weight Loss During RT') +
  xlab('Hospital') +
  ylim(c(-22.5,22.5)) +
  theme(axis.title.x = element_text(face="bold", colour="black",size=18),axis.text.x=element_text(size=20,face="bold",colour="black"))+
  theme(axis.title.y = element_text(face="bold", colour="black",size=18),axis.text.y = element_text(face="bold", colour="black",size=20))