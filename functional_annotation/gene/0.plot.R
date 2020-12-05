# Title     : TODO
# Objective : TODO
# Created by: donggy
# Created on: 2020/6/11

library(ggplot2)
library(reshape2)
library(dplyr) 



data = read.csv('a.csv')

df = data.frame(Tissues=c('1', '2', '3', '4'), # 1:1 2:2-5 3:6-42 4:43-53
                Nonassociated=as.numeric(unlist(data['Non.associated'])), 
                Only_disease=as.numeric(unlist(data['Only.disease'])),
                Comorbidity=as.numeric(unlist(data['Comorbidity'])))



# melt转换为长表格为ggplot2绘图通用格式
# geom_segment添加直线和曲线，arrange按门水平名称字母降序排列；cumsum先将数值累计，再用mutate取代；现在己有两组间的高度位置，再设置X轴位置1.25, 1.75, 和Y位置
pdf('a.pdf')
ggplot(melt(df), aes(x=variable, y=value, fill=Tissues)) + 
  geom_bar(stat = "identity", width=0.5, col='black')  + theme_classic()+
  geom_segment(data=df %>% arrange(by=desc(Tissues)) %>% mutate(GroupA=cumsum(Nonassociated)) %>% mutate(GroupB=cumsum(Only_disease)), aes(x=1.25, xend=1.75, y=GroupA, yend=GroupB))+ 
  geom_segment(data=df %>% arrange(by=desc(Tissues)) %>% mutate(GroupB=cumsum(Only_disease)) %>% mutate(GroupC=cumsum(Comorbidity)), aes(x=2.25, xend=2.75, y=GroupB, yend=GroupC))+
  labs(x=NULL,y='Gene ratio', family="Helvetica")
# 添加theme_classic()修改主题样式，这个经典主题我更喜欢
# x和xend分别为起始和终止，1，2组间X值起始分别为1.25和1.75，2，3组间则为2.25和2.75
dev.off()




data = read.csv('a1.csv')

df = data.frame(Tissues=c('1', '2', '3', '4'), # 1:1 2:2-5 3:6-42 4:43-53
                Nonassociated=as.numeric(unlist(data['Non.associated'])), 
                Only_disease=as.numeric(unlist(data['Only.disease'])),
                Comorbidity=as.numeric(unlist(data['Comorbidity'])))



# melt转换为长表格为ggplot2绘图通用格式
# geom_segment添加直线和曲线，arrange按门水平名称字母降序排列；cumsum先将数值累计，再用mutate取代；现在己有两组间的高度位置，再设置X轴位置1.25, 1.75, 和Y位置
pdf('a1.pdf')
ggplot(melt(df), aes(x=variable, y=value, fill=Tissues)) + 
  geom_bar(stat = "identity", width=0.5, col='black')  + theme_classic()+
  geom_segment(data=df %>% arrange(by=desc(Tissues)) %>% mutate(GroupA=cumsum(Nonassociated)) %>% mutate(GroupB=cumsum(Only_disease)), aes(x=1.25, xend=1.75, y=GroupA, yend=GroupB))+ 
  geom_segment(data=df %>% arrange(by=desc(Tissues)) %>% mutate(GroupB=cumsum(Only_disease)) %>% mutate(GroupC=cumsum(Comorbidity)), aes(x=2.25, xend=2.75, y=GroupB, yend=GroupC))+
  labs(x=NULL,y='Gene ratio', family="Helvetica")
# 添加theme_classic()修改主题样式，这个经典主题我更喜欢
# x和xend分别为起始和终止，1，2组间X值起始分别为1.25和1.75，2，3组间则为2.25和2.75
dev.off()