#!/usr/bin/Rscript
argv <- commandArgs(TRUE)



#Novel sequence simulation result plotting
#Duan Zhongqu, 2018-10-21


#packages 
library(ggplot2)

data<-read.table(argv[1],sep="\t",header = T)
data_dir<-argv[2]
library(ggplot2)
plot1<-paste(data_dir,"simNovelSeq.pdf",sep="")

p<-ggplot(data = data,aes(x=Sample.number,y=Sequence.length/1000000))
p=p+
  geom_point(col="red",size=1.5)+xlab("Number of individuals")+ylab("Total length (Mb)")+theme_bw()+
  theme(axis.text=element_text(size=20,face="bold",colour = "black"),axis.title=element_text(size=20,face="bold",colour = "black"))+
  xlim(0,max(data$Sample.number))+ylim(0,max(data$Sequence.length)/1000000)
ggsave(plot1,p)
