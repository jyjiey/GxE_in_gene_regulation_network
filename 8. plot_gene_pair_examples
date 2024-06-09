library(ggpubr)
matrixbd<-read.csv("~/Downloads/matrixbd_for_groupwiseWGCNA.csv", header=T,row.names=1)
RNA_Brachy_other <-read.csv("~/Desktop/codes\ jwafs/Data_otherinfo_updatedbatch.csv", header=T,row.names=1)
glucose<-log(RNA_Brachy_other$glucose.nmol.Glc.equivalents.mg)
size1=6/.pt
data <- data.frame(sample =c("Bd21d", "Bd21d", "Bd21c", "Bd21c", "Bd3-1d", "Bd3-1d","Bd3-1c", "Bd3-1c"),
                   regulator=c(1, 2, 1, 2, 1, 2,1, 2),
                   target=c(1.04,1.04,1.1,1.1,1.16,1.16,1.2,2.2))              
          ill1<- ggplot(data, aes(regulator, target, group = sample,color= sample)) + 
    geom_line(size=1)  +scale_color_manual(values=c('red', 'blue', 'orange',"purple"))  +theme(       axis.text.y=element_blank(),text = element_text(size = 6),
        axis.ticks.y=element_blank(),
       axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(color="black", size=6, face="bold")) + ggtitle("GxE slope") 

data <- data.frame(sample =c("Bd21d", "Bd21d", "Bd21c", "Bd21c", "Bd3-1d", "Bd3-1d","Bd3-1c", "Bd3-1c"),
                   regulator=c(1, 2, 1, 2, 1, 2,1, 2),
                   target=c(1.06,1.06,1.1,2.1,1.14,1.14,1.2,2.2))              
          ill2<- ggplot(data, aes(regulator, target, group = sample,color= sample)) + 
    geom_line(size=1)  +scale_color_manual(values=c('red', 'blue', 'orange',"purple"))  +theme(       axis.text.y=element_blank(),text = element_text(size = 6),
        axis.ticks.y=element_blank(),
       axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(color="black", size=6, face="bold")) + ggtitle("E slope") 


data <- data.frame(sample =c("Bd21d", "Bd21d", "Bd21c", "Bd21c", "Bd3-1d", "Bd3-1d","Bd3-1c", "Bd3-1c"),
                   regulator=c(1, 2, 1, 2, 1, 2,1, 2),
                   target=c(1.04,1.04,1.1,1.1,1.1,2.1,1.2,2.2))              
          ill3<- ggplot(data, aes(regulator, target, group = sample,color= sample)) + 
    geom_line(size=1)  +scale_color_manual(values=c('red', 'blue', 'orange',"purple"))  +theme(       axis.text.y=element_blank(),text = element_text(size = 6),
        axis.ticks.y=element_blank(),
       axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(color="black", size=6, face="bold")) + ggtitle("G slope") 

data <- data.frame(sample =c("Bd21d", "Bd21d", "Bd21c", "Bd21c", "Bd3-1d", "Bd3-1d","Bd3-1c", "Bd3-1c"),
                   regulator=c(1, 2, 1, 2, 1, 2,1, 2),
                   target=c(1.,3.,1.1,3.1,1.,1.5,1.1,1.6))              
          ill4<- ggplot(data, aes(regulator, target, group = sample,color= sample)) + 
    geom_line(size=1)  +scale_color_manual(values=c('red', 'blue', 'orange',"purple"))  +theme(       axis.text.y=element_blank(),text = element_text(size = 6),
        axis.ticks.y=element_blank(),
       axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(color="black", size=6, face="bold")) + ggtitle("G slope \n (no correlation \n change)") 

data <- data.frame(sample =c("Bd21d", "Bd21d", "Bd21c", "Bd21c", "Bd3-1d", "Bd3-1d","Bd3-1c", "Bd3-1c"),
                   regulator=c(1, 2, 1, 2, 1, 2,1, 2),
                   target=c(1.,2.,1.1,2.1,1.3,2.3,1.2,2.2))              
          ill5<- ggplot(data, aes(regulator, target, group = sample,color= sample)) + 
    geom_line(size=1)  +scale_color_manual(values=c('red', 'blue', 'orange',"purple"))  +theme(       axis.text.y=element_blank(),text = element_text(size = 6),
        axis.ticks.y=element_blank(),
       axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(color="black", size=6, face="bold")) + ggtitle("transmission") 

data <- data.frame(sample =c("Bd21d", "Bd21d", "Bd21c", "Bd21c", "Bd3-1d", "Bd3-1d","Bd3-1c", "Bd3-1c"),
                   regulator=c(1, 2, 1, 2, 1, 2,1, 2),
                   target=c(2.,3.,2.1,3.1,1.,2.,1.2,2.2))              
          ill6<- ggplot(data, aes(regulator, target, group = sample,color= sample)) + 
    geom_line(size=1)  +scale_color_manual(values=c('red', 'blue', 'orange',"purple"))  +theme(       axis.text.y=element_blank(),text = element_text(size = 6),
        axis.ticks.y=element_blank(),
       axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(color="black", size=6, face="bold"))+ ggtitle("intercept") 
ill7<-NULL 



text1 <- ggparagraph(text = "predicted relationships", face = "bold", size = 6, color = "black") +
          theme(plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1),"lines"))
text1_1 <- ggparagraph(text = "glucose gene regulation", face = "bold", size = 6, color = "black") +
          theme(plot.margin = unit(c(t = 1, r = 1, b = 1, l = 1),"lines"))
illustration_raw<-ggarrange(text1,ill1, ill2,text1_1,   ill3, ill4, ill5, ill6, labels=c("","A","B","C","D","E","F","G"),widths=c(0.6,1,1,1,1,1,1,1),ncol=8,nrow=1, common.legend=TRUE, legend="top")



names2<-"Bradi3g00920.1"
								names1<-"Bradi2g15940.1"
								
Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)==names2)])	

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
a1<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(paste0(substr(names1,1,12),"")) + ylab(paste0(substr(names2,1,12),"")) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))
condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

a2<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=paste0(substr(names2,1,12),""),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("transcript abundance")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())

names2<-"Bradi4g40850.2"
								names1<-"Bradi3g50220.1"
								
Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)==names2)])	

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
b1<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(paste0(substr(names1,1,12),"")) + ylab(paste0(substr(names2,1,12),"")) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))
condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

b2<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("transcript abundance")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())
names2<-"Bradi1g15840.1"
								names1<-"Bradi3g44990.1"
								
Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)==names2)])	

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
c1<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(paste0(substr(names1,1,12),"")) + ylab(paste0(substr(names2,1,12),"")) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))
condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

c2<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("transcript abundance")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())


names2<-"Bradi4g20520.1"
								names1<-"Bradi3g19730.4"
								
Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)==names2)])	

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
d1<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(paste0(substr(names1,1,12),"")) + ylab(paste0(substr(names2,1,12),"")) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))
condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

d2<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("transcript abundance")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())


names2<-"Bradi1g53680.1"
								names1<-"Bradi5g18277.1"								
Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)==names2)])	

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
e1<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(paste0(substr(names1,1,12),"")) + ylab(paste0(substr(names2,1,12),"")) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))
condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

e2<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("transcript abundance")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())


names2<-"Bradi3g45230.1"
								names1<-"Bradi2g10880.2"
								
								Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)== names2)])	


type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		


sample(ggplot2)
f1<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(paste0(substr(names1,1,12),"")) + ylab(paste0(substr(names2,1,12),"")) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))
condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))

gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

f2<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("transcript abundance")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())


names2<-"glucose"

								names1<-"Bradi4g39520.1"
								
								Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-log(RNA_Brachy_other$glucose.nmol.Glc.equivalents.mg)[match(rownames(matrixbd),paste0("X",RNA_Brachy_other[,1]))]

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
g1<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(paste0(substr(names1,1,12),"")) + ylab(paste0(substr(names2,1,12),"")) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))
condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

g2<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("transcript abundance \n or leaf glucose level")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())



actual_raw<-ggarrange(text_grob("actual \n cases", color = "black",face="bold", size = 6) ,a1,b1,g1,c1,d1,e1,f1,widths=c(0.6,1,1,1,1,1,1,1),ncol=8,nrow=1)



expression_raw<-ggarrange(text_grob("regulator \n and \n target \n gene \n expression", color = "black",face="bold", size = 6)  ,a2,b2,g2,c2,d2,e2,f2,widths=c(0.6,1,1,1,1,1,1,1),ncol=8,nrow=1)

data_table<-cbind(cis, cisxe, cis_c, cis_d,trans,transxe,trans_c,trans_d)
table.a<-NULL
table.b<-NULL
table.c<-NULL
table.d<-NULL
table.e<-NULL
table.f<-NULL
table.g<-NULL
table.h<-NULL
table.i<-NULL
table.j<-NULL
table.k<-NULL
table.l<-NULL
table.m<-NULL
table.n<-NULL
table.o<-NULL






	name1 <-"Bradi1g67730"
text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.b<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
	name1 <-"Bradi2g49912"
text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.c<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
                                      

                                      
	name1 <-"Bradi4g40850"
text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.e<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
	name1 <-"Bradi4g24650"
text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.f<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
                                      
                                      	name1 <-"Bradi2g52260"


	name1 <-"Bradi1g15840"
text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.h <-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6)))                          


                               
	name1 <-"Bradi3g22697"

text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.i<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 


	name1 <-"Bradi1g53680"

text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.k<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
                                            
                                     

	                                                   
	name1 <-"Bradi3g45230"

text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.l<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
                           

	                                                   
	name1 <-"Bradi3g27281"

text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.m<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
                                      
                                      
  table.m<-NULL       
                                                  	
	name1 <-"Bradi3g00920"
text.a<-rbind(data_table[which(rownames(data_table)==name1),1:4],data_table[which(rownames(data_table)==name1),5:8])
colnames(text.a)<-c("-","xE","c","d")
rownames(text.a)<-c("cis","trans")
text.a[which(text.a==0)]<-""
text.a[which(text.a==1)]<-"*"
table.a<-ggtexttable(text.a,theme = ttheme(colnames.style = colnames_style(color = "Black", fill = "white",size = 6),rownames.style = rownames_style(color = "Black", fill = "white",size = 6),
                                      tbody.style = tbody_style(color = "black", fill = "grey",size = 6))) 
ase_raw<-ggarrange(text_grob("ASE \n of \n target \n gene", color = "black",face="bold", size = 6) ,table.a, table.e,  table.o ,table.h, table.o,table.k,table.l, widths=c(0.6,1,1,1,1,1,1,1),ncol=8,nrow=1)




pdf("~/Downloads/figure6_illustrated_examples.pdf", width=10, height =6)
ggarrange(illustration_raw, actual_raw, expression_raw,ase_raw,ncol=1,nrow=4, common.legend=FALSE,align = "v",heights=c(1.5,1.5,1.5,1))



dev.off()




names1<-"Bradi3g14946.2"
								names2<-"Bradi3g33070.2"
		size1=3						
Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)==names2)])	

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
a<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(substr(names1,1,12)) + ylab(substr(names2,1,12)) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))



condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

b<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("gene expression")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())




names1<-"Bradi1g44480.1"
								names2<-"Bradi1g37410.1"
		size1=3						
Gene1<-as.numeric(matrixbd[,which(colnames(matrixbd)==names1)])
Gene2<-as.numeric(matrixbd[,which(colnames(matrixbd)==names2)])	

type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
	Data<-cbind(Gene1,Gene2,type)
	write.csv(Data,"~/Data_temperatl.csv")
	Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)		

sample(ggplot2)
c<-ggplot(Data, aes(x= Gene1, y= Gene2,
	colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48))) )+ geom_point(size = 0.8,colour=c(rep("blue",48),rep("red",48),rep("purple",48),rep("orange",48)))+ xlab(substr(names1,1,12)) + ylab(substr(names2,1,12)) + geom_smooth(data=subset(Data, type== "Bd21c" ), method='lm',colour="blue",  se=TRUE)+ geom_smooth(data=subset(Data, type== "Bd21d" ), method='lm',colour="red",  se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1c" ), method='lm',colour="purple", se= TRUE)+ geom_smooth(data=subset(Data, type== "Bd3-1d" ), method='lm', colour="orange", se=TRUE)+ theme_bw()+theme(text = element_text(size = 6))



condition <-c(rep("c",48), rep("d",48), rep("c",48), rep("d",48),rep("c",48), rep("d",48), rep("c",48), rep("d",48))
genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))
type<-c(type, type)
value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value, type)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
color<-c(c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)),c(rep("red",48),rep("blue",48),rep("orange",48),rep("purple",48)))

d<-ggplot(Data, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(value)-0.5,max(value)+2)+annotate("text",1,max(na.omit(value))+1,label=substr(names1,1,12),color="black",size=size1,hjust= 0.2,vjust=0.5)+annotate("text",1,max(na.omit(value))+1,label=substr(names2,1,12),color="grey",size= size1,hjust= 0.2,vjust=2)+ylab("gene expression")+ theme(panel.grid.major = element_blank(),text = element_text(size = 6), panel.grid.minor = element_blank(),
panel.background = element_blank())


pdf("~/Downloads/twomoreexamples.pdf", width=3, height =3)


ggarrange(c,a,d,b, labels=c("A","B","",""),ncol=2,nrow=2, common.legend=TRUE, legend="bottom",heights = c(1, 1.1,0.4,1, 1.1,0.4,1, 1.1,0.4,1, 1.1,0.4))



dev.off()

