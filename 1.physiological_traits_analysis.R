#figure1 (used full model for significant test shown in figure 1 (this code), but use the selected model for the text see supplemental table 4)
#table1 is removed from text
#figure 2
#figure 3 in supplemental data 

#section_1. Physiology plots
#require(car)
##using for distribution fit
require(MASS)
#library(MCMCglmm)
#library(scapeMCMC)
#library(agridat)
require(ggplot2)
require(lme4)
library(lmerTest)
library(MuMIn)
library("emmeans")

library(ggpubr)
##get distribution of data
theme_set(
  theme_bw()
)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

###distribution check
fitData <- function(data, fit="gamma", sample=0.5){
 distrib = list()
 numfit <- length(fit)
 results = matrix(0, ncol=5, nrow=numfit)

 for(i in 1:numfit){
if((fit[i] == "gamma") | 
     (fit[i] == "poisson") | 
     (fit[i] == "weibull") | 
     (fit[i] == "exponential") |
     (fit[i] == "logistic") |
     (fit[i] == "normal") | 
     (fit[i] == "geometric")
) 
  distrib[[i]] = fit[i]
else stop("Provide a valid distribution to fit data" )
 }

 # take a sample of dataset
 n = round(length(data)*sample)
 data = sample(data, size=n, replace=F)

 for(i in 1:numfit) {
  if(distrib[[i]] == "gamma") {
  gf_shape = "gamma"
  fd_g <- fitdistr(data, "gamma")
  est_shape = fd_g$estimate[[1]]
  est_rate = fd_g$estimate[[2]]

  ks = ks.test(data, "pgamma", shape=est_shape, rate=est_rate)

  # add to results
  results[i,] = c(gf_shape, est_shape, est_rate, ks$statistic, ks$p.value)
}

else if(distrib[[i]] == "poisson"){
  gf_shape = "poisson"
  fd_p <- fitdistr(data, "poisson")
  est_lambda = fd_p$estimate[[1]]

  ks = ks.test(data, "ppois", lambda=est_lambda)
  # add to results
  results[i,] = c(gf_shape, est_lambda, "NA", ks$statistic, ks$p.value)

}

else if(distrib[[i]] == "weibull"){
  gf_shape = "weibull"
  fd_w <- fitdistr(data,densfun=dweibull,start=list(scale=1,shape=2))
  est_shape = fd_w$estimate[[1]]
  est_scale = fd_w$estimate[[2]]

  ks = ks.test(data, "pweibull", shape=est_shape, scale=est_scale)
  # add to results
  results[i,] = c(gf_shape, est_shape, est_scale, ks$statistic, ks$p.value) 
}

else if(distrib[[i]] == "normal"){
  gf_shape = "normal"
  fd_n <- fitdistr(data, "normal")
  est_mean = fd_n$estimate[[1]]
  est_sd = fd_n$estimate[[2]]

  ks = ks.test(data, "pnorm", mean=est_mean, sd=est_sd)
  # add to results
  results[i,] = c(gf_shape, est_mean, est_sd, ks$statistic, ks$p.value)
}

else if(distrib[[i]] == "exponential"){
  gf_shape = "exponential"
  fd_e <- fitdistr(data, "exponential")
  est_rate = fd_e$estimate[[1]]
  ks = ks.test(data, "pexp", rate=est_rate)
  # add to results
  results[i,] = c(gf_shape, est_rate, "NA", ks$statistic, ks$p.value)
}

else if(distrib[[i]] == "logistic"){
  gf_shape = "logistic"
  fd_l <- fitdistr(data, "logistic")
  est_location = fd_l$estimate[[1]]
  est_scale = fd_l$estimate[[2]]
  ks = ks.test(data, "plogis", location=est_location, scale=est_scale)
  # add to results
  results[i,] = c(gf_shape, est_location, est_scale, ks$statistic,    ks$p.value) 
    }
  }
#this function is modified from online resources
  results = rbind(c("distribution", "param1", "param2", "ks stat", "ks    pvalue"),   results)
  #print(results)
  return(results)
  }
head(Drydown_brachy)  

head(Drydown_brachy)
matrix<-matrix(0,17,7)
colnames(matrix)<-c("GxExT","GxE","GxT","ExT","T","E","G")
rownames(matrix)[1]<-"A"
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
m=1

  pdf(paste0("~/Downloads/figure1_dry-down.pdf"), width=14, height = 8)

###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-Drydown_brachy $water.usage
#3.Netphotosynthesisrate3$A<-c(-0.1*Drydown_brachy $hydraulic.potential.1)
#4.Netphotosynthesisrate3$A<-Drydown_brachy $A
#5.Netphotosynthesisrate3$A<-Drydown_brachy $gsw
#6.Netphotosynthesisrate3$A<-(Drydown_brachy $fresh.shoot+Drydown_brachy $fresh.root)
#7.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment | Accession *harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),5), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale



g1<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Water usage",~'(g/day)' ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1, group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2)+facet_grid(rows=vars(Accession))


###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-c(-0.1*Drydown_brachy $hydraulic.potential.1)
#4.Netphotosynthesisrate3$A<-Drydown_brachy $A
#5.Netphotosynthesisrate3$A<-Drydown_brachy $gsw
#6.Netphotosynthesisrate3$A<-(Drydown_brachy $fresh.shoot+Drydown_brachy $fresh.root)
#7.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment | Accession *harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale


g2<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Hydraulic potential",~'(MPa)' ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))
  
 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-Drydown_brachy $A
#5.Netphotosynthesisrate3$A<-Drydown_brachy $gsw
#6.Netphotosynthesisrate3$A<-(Drydown_brachy $fresh.shoot+Drydown_brachy $fresh.root)
#7.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment | Accession *harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale

 
  
  
##rerun
g3<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Net Carbon Assimilation",~'(',mu,'(mol',~CO[2],~m^{-2},~s^{-1},')' ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1, group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))


 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-Drydown_brachy $gsw
#6.Netphotosynthesisrate3$A<-(Drydown_brachy $fresh.shoot+Drydown_brachy $fresh.root)
#7.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment |Accession*harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale

 

g4<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("gsw",~'(mol',~CO[2],~m^{-2},~s^{-1},')' ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1, group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))

##rerun
 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-(Drydown_brachy $fresh.shoot+Drydown_brachy $fresh.root)
#7.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment | Accession *harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale

 

g5<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("fresh weight",~'(g)' ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))

##rerun
 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment | Accession *harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale

 



g6<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("dry weight",~'(g)' ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))

 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment | Accession *harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale

 

g7<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("dry root shoot ratio" ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))




setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "glucose.nmol.Glc.equivalents.mg"))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$glucose.nmol.Glc.equivalents.mg, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $glucose.nmol.Glc.equivalents.mg)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- glm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3,family=Gamma(link = "inverse"))
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment +Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
b<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment + Treatment*harvest.day+(1| freeze.dry.batch), data= glucose3,family=Gamma(link = "inverse"))
c<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment + Treatment*harvest.day+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
d<- glm(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment +Treatment*harvest.day, data= glucose3,family=Gamma(link = "inverse"))
anova(a,b,c,d)

glucose3 $harvest.day<-as.factor(glucose3 $harvest.day)
d<- glmer(glucose ~  Accession* Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
gm_mc <- emmeans(d, ~ Treatment |Accession* harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 

letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)


 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g8<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Glucose",~'(','nmol /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))

###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
head(Drydown_brachy)

Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-Drydown_brachy $NPQ
Netphotosynthesisrate3$A[which(Netphotosynthesisrate3$A<0)]<-NA
#3.Netphotosynthesisrate3$A<-c(-0.1*Drydown_brachy $hydraulic.potential.1)
#4.Netphotosynthesisrate3$A<-Drydown_brachy $A
#5.Netphotosynthesisrate3$A<-Drydown_brachy $gsw
#6.Netphotosynthesisrate3$A<-(Drydown_brachy $fresh.shoot+Drydown_brachy $fresh.root)
#7.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))
#anova(c1)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
c<- lm(Netphotosynthesisrate ~ Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
summary(c)
gm_mc <- emmeans(c, ~ Treatment |Accession* harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 
letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale



g9<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("NPQ" ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1, group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2)+facet_grid(rows=vars(Accession))

###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
head(Drydown_brachy)

Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]
Netphotosynthesisrate3$A<-Drydown_brachy $Average.of.Relative.Chlorophyll
Netphotosynthesisrate3$A[which(Netphotosynthesisrate3$A<0)]<-NA
#3.Netphotosynthesisrate3$A<-c(-0.1*Drydown_brachy $hydraulic.potential.1)
#4.Netphotosynthesisrate3$A<-Drydown_brachy $A
#5.Netphotosynthesisrate3$A<-Drydown_brachy $gsw
#6.Netphotosynthesisrate3$A<-(Drydown_brachy $fresh.shoot+Drydown_brachy $fresh.root)
#7.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.shoot+Drydown_brachy $dry.root)
#8.Netphotosynthesisrate3$A<-(Drydown_brachy $dry.root/Drydown_brachy $dry.shoot)
#9.Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
res = fitData(na.omit(Netphotosynthesisrate3$A), fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)

res
###
Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)
options(na.action = "na.fail")
fm1 <- lm(A~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = Netphotosynthesisrate3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1



#c2<- lm(Netphotosynthesisrate ~   harvest.day+ Treatment, data= Netphotosynthesisrate3)
#c1<- lm(Netphotosynthesisrate ~  Accession*Treatment*harvest.day, data= Netphotosynthesisrate3)
#anova(c2,c1)
#summary(anova(c1,c2))

letter <- rep( " ",6)        
theme_set(
  theme_bw()
)

letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right

scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale




g10<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Relative chlorophyll content" ))))+xlab("Time (day)") + geom_text(
    aes(label = letterAAA, y= letterlocation1, group=Treatment),color="black", position = position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2)+facet_grid(rows=vars(Accession))

  pdf(paste0("~/Downloads/figure1_dry-down.pdf"), width=10, height = 6)

  #pdf(paste0("~/Downloads/Photosynthesis_normalized.pdf"), width=8, height = 5)
# ggarrange(g1, g2,g3,g4,g8,labels=c("A","B","C","D","E"),ncol=5,nrow=1, common.legend=TRUE, legend="bottom")
 ggarrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,labels=c("a","b","c","d","e","f","g","h","i","j"),ncol=5,nrow=2, common.legend=TRUE, legend="bottom")
#ggarrange(g1,g9,g3,g4,g8,g2,g5,g6,g7,g10,labels=c("A","B","C","D","E","F","G","H","I","J"),ncol=5,nrow=2, common.legend=TRUE, legend="bottom")
dev.off()


letter <- symnum(matrix, corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        

textMatrix = matrix(paste(  signif(matrix, 2),letter),17,7)
write.csv(textMatrix,"~/Downloads/table1.csv")


pd <- position_dodge(0.8) # move them .05 to the left and right
RNA_Brachy_other <-read.csv("~/Desktop/codes jwafs/Data_otherinfo_updatedbatch.csv", header=T)
head(RNA_Brachy_other)
RNA_Brachy_other<-RNA_Brachy_other[-c(34,210),]
RNA_Brachy_other$AccessionTreatment<-paste0(RNA_Brachy_other$Accession, RNA_Brachy_other$Treatment)
RNA_Brachy_other$AccessionTreatment<-as.factor(RNA_Brachy_other$AccessionTreatment)
a<-RNA_Brachy_other$ day5_waterusage
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 
res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment+(1| batch1experiment), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment, data= data)
anova(a,b)
a<-a
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
unit<-expression(plain(paste("Water Usage",~'(g)' )))

#annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(anova(a)[1,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[1,6]<2.2e-16, 2.2e-16,anova(a)[1,6]),2)),label2=paste0("italic(P(E))",ifelse(anova(a)[2,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[2,6]<2.2e-16, 2.2e-16,anova(a)[2,6]),2)),label3=paste0("italic(P(GxE))",ifelse(anova(a)[3,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[3,6]<2.2e-16, 2.2e-16,anova(a)[3,6]),2)))

full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e1<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,1.3),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(unit)  + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


head(RNA_Brachy_other)
a<-RNA_Brachy_other$ water_content*100
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment+(1| batch1experiment), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment, data= data)
anova(a,b)
a<-a
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.5* scale
unit<-expression(plain(paste("Water Content" ,~'(%)')))

full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate)+0.3* scale,min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
#annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(anova(a)[1,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[1,6]<2.2e-16, 2.2e-16,anova(a)[1,6]),2)),label2=paste0("italic(P(E))",ifelse(anova(a)[2,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[2,6]<2.2e-16, 2.2e-16,anova(a)[2,6]),2)),label3=paste0("italic(P(GxE))",ifelse(anova(a)[3,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[3,6]<2.2e-16, 2.2e-16,anova(a)[3,6]),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e2<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/50, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.45,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(unit)  + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


a<-RNA_Brachy_other$ Photosynthesis_A
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment+(1| batch1experiment), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment, data= data)
anova(a,b)
a<-a
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale

full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))

#annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(anova(a)[1,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[1,6]<2.2e-16, 2.2e-16,anova(a)[1,6]),2)),label2=paste0("italic(P(E))",ifelse(anova(a)[2,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[2,6]<2.2e-16, 2.2e-16,anova(a)[2,6]),2)),label3=paste0("italic(P(GxE))",ifelse(anova(a)[3,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[3,6]<2.2e-16, 2.2e-16,anova(a)[3,6]),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e3<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,1.3),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Net Carbon Assimilation",~'(',mu,'mol',~CO[2],~m^{-2},~s^{-1},')' ))))  + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)

head(RNA_Brachy_other)
a<-RNA_Brachy_other$ Photosynthesis_gsw
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment+(1| batch1experiment), data= data,family=Gamma(link = "inverse"))
b<- glm(Netphotosynthesisrate ~ Accession+ Treatment+Accession* Treatment, data= data,family=Gamma(link = "inverse"))
anova(a,b)
a<-a

##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale

full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,min(data $Netphotosynthesisrate),max(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
#Anova(a,type="III")
#annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(anova(a)[1,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[1,6]<2.2e-16, 2.2e-16,anova(a)[1,6]),2)),label2=paste0("italic(P(E))",ifelse(anova(a)[2,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[2,6]<2.2e-16, 2.2e-16,anova(a)[2,6]),2)),label3=paste0("italic(P(GxE))",ifelse(anova(a)[3,6]<2.2e-16, "<","=="), signif(ifelse(anova(a)[3,6]<2.2e-16, 2.2e-16,anova(a)[3,6]),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e4<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/60, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("gsw",~'(mol',~CO[2],~m^{-2},~s^{-1},')' )))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)



head(RNA_Brachy_other)
a<-RNA_Brachy_other$ glucose.nmol.Glc.equivalents.mg
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data,family=Gamma(link = "inverse"))
b<- glm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data,family=Gamma(link = "inverse"))
c<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data,family=Gamma(link = "inverse"))
d<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data,family=Gamma(link = "inverse"))
e<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data,family=Gamma(link = "inverse"))
f<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data,family=Gamma(link = "inverse"))
g<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data,family=Gamma(link = "inverse"))
h<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data,family=Gamma(link = "inverse"))
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e5<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Glucose",~'(','nmol /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


head(RNA_Brachy_other)
a<-RNA_Brachy_other$ fructose.nmol.Glc.equivalents.mg
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data,family=Gamma(link = "inverse"))
b<- glm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data,family=Gamma(link = "inverse"))
c<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data,family=Gamma(link = "inverse"))
d<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data,family=Gamma(link = "inverse"))
e<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data,family=Gamma(link = "inverse"))
f<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data,family=Gamma(link = "inverse"))
g<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data,family=Gamma(link = "inverse"))
h<- glmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data,family=Gamma(link = "inverse"))
anova(a,b,c,d,e,f,g,h)
a<-f
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e6<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.5,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Fructose",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


head(RNA_Brachy_other)
a<-RNA_Brachy_other$ sucrose.nmol.Glc.equivalents.mg
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-e
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.4* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate)+0.3* scale,min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e7<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,2),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Sucrose",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)

head(RNA_Brachy_other)
a<-RNA_Brachy_other$ Low.DP.fructan.nmol.Glc.equivalents.mg
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 
res = fitData(data$Netphotosynthesisrate, fit=c("gamma","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.4* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate)+0.2*scale,min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e8<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,1.5),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Low.DP.fructan",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)

head(RNA_Brachy_other)
a<-RNA_Brachy_other$ Starch.nmol.Glc.equivalents.mg
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 
res = fitData(data$Netphotosynthesisrate, fit=c("gamma","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch3pellete), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch3pellete), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch3pellete), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch3pellete), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-g
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.4* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate)+0.2* scale,min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e9<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,1.5),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Starch",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


head(RNA_Brachy_other)
a<-RNA_Brachy_other$ High.DP.Fructan.nmol.Glc.equivalents.mg
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 
res = fitData(data$Netphotosynthesisrate, fit=c("gamma","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch3pellete), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch3pellete), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch3pellete), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch3pellete), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-g
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.4* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=2,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate)+0.2* scale,min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e10<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,1.5),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("High.DP.fructan",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)

head(RNA_Brachy_other)
a<-RNA_Brachy_other$ Protein...ug.mg.
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 
res = fitData(data$Netphotosynthesisrate, fit=c("gamma","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch3pellete), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch3pellete), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch3pellete), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch3pellete), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-f
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.4* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,min(data $Netphotosynthesisrate),max(data $Netphotosynthesisrate))+0.3*scale,label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e11<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,1.5),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Protein",~'(','nmol /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


head(RNA_Brachy_other)
a<-RNA_Brachy_other$ amino.acids..nmol.mg.
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other[,3],RNA_Brachy_other[,4])=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other[,3])[-outlier]
Treatment<-as.factor(RNA_Brachy_other[,4])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other[,5])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other[,20])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other[,22])[-outlier]
batch4sn<-as.factor(RNA_Brachy_other[,26])[-outlier]

data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch3pellete, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 
res = fitData(data$Netphotosynthesisrate, fit=c("gamma","normal","exponential","poisson","exponential"),
    sample=1)
res
a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,min(data $Netphotosynthesisrate),max(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))

e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e12<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =min(scale/10,1.5),
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Amino Acids",~'(','nmol /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)



pdf(paste0("~/Downloads/figure2_more.pdf"), width=9, height = 6)

ggarrange(e1,e3,e5,e6,e7,e8,e9,e10,e11,e12,labels=c("A","B","C","D","E","F","G","H","I","J"),ncol=5,nrow=2, common.legend=TRUE, legend="bottom")
dev.off()





head(Drydown_brachy)
matrix<-matrix(0,17,7)
colnames(matrix)<-c("GxExT","GxE","GxT","ExT","T","E","G")
rownames(matrix)[1]<-"A"
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
m=1
###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)



setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "glucose.nmol.Glc.equivalents.mg"))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$glucose.nmol.Glc.equivalents.mg, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $glucose.nmol.Glc.equivalents.mg)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- glm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3,family=Gamma(link = "inverse"))
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment +Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
b<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment + Treatment*harvest.day+(1| freeze.dry.batch), data= glucose3,family=Gamma(link = "inverse"))
c<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment + Treatment*harvest.day+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
d<- glm(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment +Treatment*harvest.day, data= glucose3,family=Gamma(link = "inverse"))
anova(a,b,c,d)

glucose3 $harvest.day<-as.factor(glucose3 $harvest.day)
d<- glm(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment +Treatment*harvest.day, data= glucose3,family=Gamma(link = "inverse"))
gm_mc <- emmeans(d, ~ Treatment |Accession* harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 

letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)


 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $glucose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g1<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Glucose",~'(','nmol /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))




glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "fructose.nmol.Glc.equivalents.mg"))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$fructose.nmol.Glc.equivalents.mg, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $fructose.nmol.Glc.equivalents.mg)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- glm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3,family=Gamma(link = "inverse"))
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day +Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
b<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day + Treatment*harvest.day+(1| freeze.dry.batch), data= glucose3,family=Gamma(link = "inverse"))
c<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day + Treatment*harvest.day+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
d<- glm(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day +Treatment*harvest.day, data= glucose3,family=Gamma(link = "inverse"))
anova(a,b,c,d)

glucose3 $harvest.day<-as.factor(glucose3 $harvest.day)
d<- glmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day + Treatment*harvest.day+(1| DW.Plate..), data= glucose3,family=Gamma(link = "inverse"))
gm_mc <- emmeans(d, ~ Treatment |Accession* harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 

letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)


 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $fructose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g2<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Fructose",~'(','nmol Glc equiv /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))


glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "sucrose.nmol.Glc.equivalents.mg"))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$sucrose.nmol.Glc.equivalents.mg, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $sucrose.nmol.Glc.equivalents.mg)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- lm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- lmer(glucose ~  Accession +(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3)
b<- lmer(glucose ~  Accession +(1| freeze.dry.batch), data= glucose3)
c<- lmer(glucose ~  Accession +(1| DW.Plate..), data= glucose3)
d<- lm(glucose ~  Accession , data= glucose3)
anova(a,b,c,d)
letter<-rep("",12)
 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $sucrose.nmol.Glc.equivalents.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g3<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Sucrose",~'(','nmol Glc equiv /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))





setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "Low.DP.fructan.nmol.Glc.equivalents.mg..plate2.has.batch.effect.."))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$Low.DP.fructan.nmol.Glc.equivalents.mg..plate2.has.batch.effect.., fit=c("gamma","normal"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $Low.DP.fructan.nmol.Glc.equivalents.mg..plate2.has.batch.effect..)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- lm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- lmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment +Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3)
b<- lmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment + Treatment*harvest.day+(1| freeze.dry.batch), data= glucose3)
c<- lmer(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment + Treatment*harvest.day+(1| DW.Plate..), data= glucose3)
d<- lm(glucose ~  Accession + Treatment +harvest.day+ Accession* Treatment +Treatment*harvest.day, data= glucose3)
anova(a,b,c,d)

letter<-rep("",12)

 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $Low.DP.fructan.nmol.Glc.equivalents.mg..plate2.has.batch.effect..


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g4<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("low DP fructan",~'(','nmol Glc equiv /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))




setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "starch.nmol.Glc.mg"))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$starch.nmol.Glc.mg, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $starch.nmol.Glc.mg)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- lm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- lmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day +Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3)
b<- lmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day + Treatment*harvest.day+(1| freeze.dry.batch), data= glucose3)
c<- lmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day + Treatment*harvest.day+(1| DW.Plate..), data= glucose3)
d<- lm(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day +Treatment*harvest.day, data= glucose3)
anova(a,b,c,d)

glucose3 $harvest.day<-as.factor(glucose3 $harvest.day)
d<- lmer(glucose ~  Accession + Treatment +harvest.day+ Accession* harvest.day + Treatment*harvest.day+(1| DW.Plate..), data= glucose3)
gm_mc <- emmeans(d, ~ Treatment |Accession* harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 

letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)


 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $starch.nmol.Glc.mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g5<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("Starch",~'(','nmol Glc equiv /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))



setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "High.DP.Fructan.nmol.Glc.mg.nmol.Glc.equivalents..mg"))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$High.DP.Fructan.nmol.Glc.mg.nmol.Glc.equivalents..mg, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $High.DP.Fructan.nmol.Glc.mg.nmol.Glc.equivalents..mg)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- lm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1

letter=rep("",12)
 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $High.DP.Fructan.nmol.Glc.mg.nmol.Glc.equivalents..mg


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g6<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("high DP fructan",~'(','nmol Glc equiv /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))



setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "Protein...ug.mg."))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$Protein...ug.mg., fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $Protein...ug.mg.)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- lm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- lmer(glucose ~  Accession + Treatment +harvest.day+Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3)
b<- lmer(glucose ~  Accession + Treatment +harvest.day+ Treatment*harvest.day+(1| freeze.dry.batch), data= glucose3)
c<- lmer(glucose ~  Accession + Treatment +harvest.day + Treatment*harvest.day+(1| DW.Plate..), data= glucose3)
d<- lm(glucose ~  Accession + Treatment +harvest.day+Treatment*harvest.day, data= glucose3)
anova(a,b,c,d)

glucose3 $harvest.day<-as.factor(glucose3 $harvest.day)
d<- lm(glucose ~  Accession + Treatment +harvest.day+Treatment*harvest.day, data= glucose3)
gm_mc <- emmeans(d, ~ Treatment |Accession* harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 

letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)


 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $Protein...ug.mg.


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g7<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("protein",~'(','nmol /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))





setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)
glucose3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day","freeze.dry.batch","DW.Plate..", "amino.acids..nmol.mg."))]
glucose3 <-na.omit(glucose3)
head(Drydown_brachy)
res = fitData(glucose3$amino.acids..nmol.mg., fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res
#
glucose3 $glucose<-as.numeric(glucose3 $amino.acids..nmol.mg.)
glucose3 $harvest.day<-as.numeric(glucose3 $harvest.day)
glucose3 $Accession <-as.factor(glucose3 $Accession)
glucose3 $Treatment <-as.factor(glucose3 $usage)
glucose3 $freeze.dry.batch<-as.factor(glucose3 $freeze.dry.batch)
glucose3 $DW.Plate..<-as.factor(glucose3 $DW.Plate..)
options(na.action = "na.fail")


fm1 <- lm(glucose~ Accession + Treatment +harvest.day+Accession* Treatment*harvest.day,data = glucose3)
dredge(fm1, rank="AIC",extra = c("R^2", adjRsq=function(x) summary(x)$adj.r.squared))
full_model<-fm1
no_3inter_model<-update(full_model,.~.-Accession : 
    Treatment : harvest.day)
no_gxeinter_model<-update(no_3inter_model,.~.-Accession : 
    Treatment)
no_gxtinter_model<-update(no_3inter_model,.~.-Accession : 
    harvest.day)
no_txeinter_model<-update(no_3inter_model,.~.-harvest.day : 
    Treatment)
g1_model <-update(no_3inter_model,.~.-harvest.day : 
    Accession-Accession: Treatment)
g2_model <-update(g1_model,.~.-Accession)

e1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-Accession: Treatment)
e2_model <-update(e1_model,.~.-Treatment)

t1_model <-update(no_3inter_model,.~.-harvest.day : 
    Treatment-harvest.day : 
    Accession)
t2_model <-update(t1_model,.~.-harvest.day)



matrix[m,]<-c(anova(full_model, no_3inter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_gxtinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(no_3inter_model, no_txeinter_model,test="Chisq")[2,]$"Pr(>Chi)",anova(t1_model, t2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(e1_model, e2_model,test="Chisq")[2,]$"Pr(>Chi)",anova(g1_model, g2_model,test="Chisq")[2,]$"Pr(>Chi)")
m=m+1


a<- lmer(glucose ~  Treatment +harvest.day+Treatment*harvest.day+(1| freeze.dry.batch)+(1| DW.Plate..), data= glucose3)
b<- lmer(glucose ~  Treatment +harvest.day+ Treatment*harvest.day+(1| freeze.dry.batch), data= glucose3)
c<- lmer(glucose ~  Treatment +harvest.day + Treatment*harvest.day+(1| DW.Plate..), data= glucose3)
d<- lm(glucose ~  Treatment +harvest.day+Treatment*harvest.day, data= glucose3)
anova(a,b,c,d)

glucose3 $harvest.day<-as.factor(glucose3 $harvest.day)
d<- lm(glucose ~  Treatment +harvest.day+Treatment*harvest.day, data= glucose3)
gm_mc <- emmeans(d, ~ Treatment |Accession* harvest.day , ddf="kenward-roger")
Hf_Aarea_lsmeans <-pairs(gm_mc)
Hf_Aarea_lsmeans 

letter<-Hf_Aarea_lsmeans
letter<-as.data.frame(letter)
letter <- symnum(letter[,dim(letter)[2]], corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        
theme_set(
  theme_bw()
)


 ###dry-down data
setwd("/Users/jiey/Desktop/codes jwafs")
Drydown_brachy <-read.csv("dry-down process data_bigmatrix.csv", header=T)


Netphotosynthesisrate3<-Drydown_brachy[,which(colnames(Drydown_brachy)%in%c("Accession", "usage", "harvest.day", "A"))]

Netphotosynthesisrate3$A<-Drydown_brachy $amino.acids..nmol.mg.


Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)

Netphotosynthesisrate3 $Netphotosynthesisrate<-as.numeric(na.omit(Netphotosynthesisrate3 $A))
Netphotosynthesisrate3 $harvest.day<-as.numeric(Netphotosynthesisrate3 $harvest.day)
Netphotosynthesisrate3 $Accession <-as.factor(Netphotosynthesisrate3 $Accession)
Netphotosynthesisrate3 $Treatment <-as.factor(Netphotosynthesisrate3 $usage)

Netphotosynthesisrate3 <-na.omit(Netphotosynthesisrate3)
Netphotosynthesisrate3 $harvest.day<-as.factor(Netphotosynthesisrate3 $harvest.day)
letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), letter)
#letterA<-cbind(c(0,0,2,2,4,4,5,5,6,6,7,7),rep(c("bd21","bd3-1"),6), rep(letter,each=2))

BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD21",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd21"
letterAA<-letterA[which(letterA[,2]=="bd21"),]
dim(BD21)
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc1<-tgc
BD21<-Netphotosynthesisrate3[Netphotosynthesisrate3 $Accession=="BD3-1",]
BD21 $Harvest.day <- as.factor(BD21 $harvest.day)
BD21 $Treatment<- as.factor(BD21 $Treatment)
BD21 $Accession<- "Bd3-1"
letterAA<-letterA[which(letterA[,2]=="bd3-1"),]
tgc <- summarySE(BD21, measurevar="Netphotosynthesisrate", groupvars=c("Accession","Treatment","harvest.day"))
letterAAA<-letterAA[match(tgc[,3], letterAA[,1]),3]
tgc<-cbind(tgc, letterAAA)
tgc<-rbind(tgc, tgc[1,])
tgc[12,2]<-"Drought"
tgc<-rbind(tgc1,tgc)

tgc <-cbind(tgc,paste0(tgc $Accession, "-d", tgc $harvest.day))
colnames(tgc)[dim(tgc)[2]]<-"g_day"
letterlocation<-aggregate(tgc$Netphotosynthesisrate,list(tgc$g_day),mean)
letterlocation1<-letterlocation[match( tgc $g_day, letterlocation[,1]),2]
tgc<-cbind(tgc, letterlocation1)
pd <- position_dodge(0.1) # move them .05 to the left and right
scale<-max(tgc$Netphotosynthesisrate)-min(tgc$Netphotosynthesisrate)
yscale1<-min(tgc$Netphotosynthesisrate)-0.2* scale
yscale2<-max(tgc$Netphotosynthesisrate)+0.2* scale
 
 
g8<-ggplot(tgc, aes(x= harvest.day, y= Netphotosynthesisrate, colour= Treatment)) + 
    geom_errorbar(aes(ymin= Netphotosynthesisrate-se, ymax= Netphotosynthesisrate +se), width=.1, position=pd) + geom_line(aes(group= Treatment),position=pd) +
    geom_point(position=pd)+
  scale_color_manual(values = c( "#00AFBB","#E7B800"))+ylab(expression(plain(paste("amino acid",~'(','nmol /',~'mg DW',')'))))+xlab("Time (day)")  + geom_text(
    aes(label = letterAAA, y= letterlocation1,  group=Treatment),color="black", position =  position_dodge(width=0), size=6/.pt)+ theme(text = element_text(size = 8),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+
  ylim(yscale1, yscale2) +facet_grid(rows=vars(Accession))
 pdf(paste0("~/Downloads/supplemental_figure1.pdf"), width=14, height = 8)

  #pdf(paste0("~/Downloads/Photosynthesis_normalized.pdf"), width=14, height = 8)
 ggarrange(g1, g2,g3,g4,g5,g6,g7,g8,labels=c("A","B","C","D","E","F","G","H"),ncol=4,nrow=2, common.legend=TRUE, legend="bottom")
 
#ggarrange(g1,g9,g3,g4,g8,g2,g5,g6,g7,g10,labels=c("a","b","c","d","e","f","g","h","i","j"),ncol=5,nrow=2, common.legend=TRUE, legend="bottom")
dev.off()
letter <- symnum(matrix, corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))        

textMatrix = matrix(paste(  signif(matrix, 2),letter),17,7)
write.csv(textMatrix,"~/Downloads/supplemental_table_1.csv")


RNA_Brachy_other <-read.csv("~/Desktop/codes jwafs/Data_otherinfo_updatedbatch.csv", header=T)
head(RNA_Brachy_other)
RNA_Brachy_other<-RNA_Brachy_other[-c(34,210),]
RNA_Brachy_other2 <-read.csv("~/Desktop/draft_SI/draft_data/RNAplants_stem_metabolisms_BIGMATRIX.csv", header=T)


RNA_Brachy_other2[,13]<-RNA_Brachy_other[match(RNA_Brachy_other2[,2], RNA_Brachy_other[,2]),]$Treatment
RNA_Brachy_other2[,14]<-RNA_Brachy_other[match(RNA_Brachy_other2[,2], RNA_Brachy_other[,2]),]$Accession
RNA_Brachy_other2[,15]<-RNA_Brachy_other[match(RNA_Brachy_other2[,2], RNA_Brachy_other[,2]),]$Experiment.Batch
RNA_Brachy_other2<-RNA_Brachy_other2[-c(8,31),]
colnames(RNA_Brachy_other2)[13]<-"Treatment"
colnames(RNA_Brachy_other2)[14]<-"Accession"
colnames(RNA_Brachy_other2)[15]<-"Experiment.Batch"

a<-RNA_Brachy_other2 $ glucose
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))

Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)
Accession<-as.factor(RNA_Brachy_other2[,14])
Treatment<-as.factor( RNA_Brachy_other2[,13])
batch1experiment<-as.factor(RNA_Brachy_other2[,15])
batch2dry <-as.factor(RNA_Brachy_other2[,12])
batch4sn <-as.factor(RNA_Brachy_other2[,1])



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-c
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e1<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Glucose",~'(','nmol /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)




a<-RNA_Brachy_other2 $ fructose
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))

Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)
Accession<-as.factor(RNA_Brachy_other2[,14])
Treatment<-as.factor( RNA_Brachy_other2[,13])
batch1experiment<-as.factor(RNA_Brachy_other2[,15])
batch2dry <-as.factor(RNA_Brachy_other2[,12])
batch4sn <-as.factor(RNA_Brachy_other2[,1])



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-d
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e2<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("Fructose",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


head(RNA_Brachy_other2)
a<-RNA_Brachy_other2 $ sucrose
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))

Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)
Accession<-as.factor(RNA_Brachy_other2[,14])
Treatment<-as.factor( RNA_Brachy_other2[,13])
batch1experiment<-as.factor(RNA_Brachy_other2[,15])
batch2dry <-as.factor(RNA_Brachy_other2[,12])
batch4sn <-as.factor(RNA_Brachy_other2[,1])



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e3<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("sucrose",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)



head(RNA_Brachy_other2)
a<-RNA_Brachy_other2 $ low.DP.fructan
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
outlier
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)[-outlier]
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch4sn <-as.factor(RNA_Brachy_other2[,1])[-outlier]



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e4<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("low DP fructan",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)


head(RNA_Brachy_other2)
a<-RNA_Brachy_other2 $ starch.nmol.Glc.mg
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
outlier
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)
Accession<-as.factor(RNA_Brachy_other2[,14])
Treatment<-as.factor( RNA_Brachy_other2[,13])
batch1experiment<-as.factor(RNA_Brachy_other2[,15])
batch2dry <-as.factor(RNA_Brachy_other2[,12])
batch4sn <-as.factor(RNA_Brachy_other2[,1])



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-c
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e5<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("starch",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)

head(RNA_Brachy_other2)
a<-RNA_Brachy_other2 $ high.DP.Fructan
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
outlier
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)
Accession<-as.factor(RNA_Brachy_other2[,14])
Treatment<-as.factor( RNA_Brachy_other2[,13])
batch1experiment<-as.factor(RNA_Brachy_other2[,15])
batch2dry <-as.factor(RNA_Brachy_other2[,12])
batch4sn <-as.factor(RNA_Brachy_other2[,1])



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e6<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("high DP fructan",~'(','nmol Glc equiv /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)



head(RNA_Brachy_other2)
a<-RNA_Brachy_other2 $ protein...ug.mg.
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
outlier
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)
Accession<-as.factor(RNA_Brachy_other2[,14])
Treatment<-as.factor( RNA_Brachy_other2[,13])
batch1experiment<-as.factor(RNA_Brachy_other2[,15])
batch2dry <-as.factor(RNA_Brachy_other2[,12])
batch4sn <-as.factor(RNA_Brachy_other2[,1])



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e7<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("protein",~'(','nmol /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)



head(RNA_Brachy_other2)
a<-RNA_Brachy_other2 $ amino.acids
###fructose find outliers,RNA anova, and boxplot
outlier<-c(match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 normal")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD21 drought")], plot=FALSE)$out,a),match(boxplot(a[which(paste(RNA_Brachy_other2 $Accession, RNA_Brachy_other2 $Treatment)=="BD3-1 drought")], plot=FALSE)$out,a),which(is.na(a)))
outlier
Netphotosynthesisrate<-as.vector(a[-outlier])
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch3pellete <-as.factor(RNA_Brachy_other2[,1])[-outlier]

Netphotosynthesisrate<-as.vector(a)[-outlier]
Accession<-as.factor(RNA_Brachy_other2[,14])[-outlier]
Treatment<-as.factor( RNA_Brachy_other2[,13])[-outlier]
batch1experiment<-as.factor(RNA_Brachy_other2[,15])[-outlier]
batch2dry <-as.factor(RNA_Brachy_other2[,12])[-outlier]
batch4sn <-as.factor(RNA_Brachy_other2[,1])[-outlier]



data<-as.data.frame(cbind(Netphotosynthesisrate, Accession, Treatment, batch1experiment, batch2dry, batch4sn))
data[data$Accession==1,2]<-rep("Bd21", length(data[data$Accession==1,1]))
data[data$Accession==2,2]<-rep("Bd3-1", length(data[data$Accession==2,1]))
data[data$Treatment==1,3]<-rep("Drought", length(data[data$Treatment ==1,1]))
data[data$Treatment ==2,3]<-rep("Control", length(data[data$Treatment ==2,1]))
#distribution 

res = fitData(data$Netphotosynthesisrate, fit=c("gamma","logistic","normal","exponential","poisson","exponential"),
    sample=1)
res


a<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry)+(1| batch4sn), data= data)
b<- lm(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession, data= data)
c<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment), data= data)
d<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession+(1| batch2dry), data= data)
e<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch4sn), data= data)
f<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch2dry)+(1| batch4sn), data= data)
g<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch4sn), data= data)
h<- lmer(Netphotosynthesisrate ~ Accession+ Treatment+ Treatment*Accession +(1| batch1experiment)+(1| batch2dry), data= data)
anova(a,b,c,d,e,f,g,h)
a<-b
##plot
scale<-max(data$Netphotosynthesisrate)-min(data$Netphotosynthesisrate)
yscale1<-min(data $Netphotosynthesisrate)-0.2* scale
yscale2<-max(data $Netphotosynthesisrate)+0.2* scale
full_model<-a
no_inter_model<-update(full_model,.~.-Accession: Treatment)
G_term_model<-update(no_inter_model,.~.-Accession)
E_term_model<-update(no_inter_model,.~.-Treatment)
interaction_pVALUE<-c(anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(full_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
G_pVALUE<-c(anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(G_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
E_pVALUE<-c(anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chisq)",anova(E_term_model, no_inter_model,test="Chisq")[2,]$"Pr(>Chi)")
annotation<-data.frame(x=1,y=ifelse(summary(a)$coefficients[3,1]<0,max(data $Netphotosynthesisrate),min(data $Netphotosynthesisrate)),label1=paste0("italic(P(G))",ifelse(G_pVALUE <2.2e-16, "<","=="), signif(ifelse(G_pVALUE <2.2e-16, 2.2e-16, G_pVALUE),2)),label2=paste0("italic(P(E))",ifelse(E_pVALUE <2.2e-16, "<","=="), signif(ifelse(E_pVALUE <2.2e-16, 2.2e-16, E_pVALUE),2)),label3=paste0("italic(P(GxE))",ifelse(interaction_pVALUE <2.2e-16, "<","=="), signif(ifelse(interaction_pVALUE <2.2e-16, 2.2e-16, interaction_pVALUE),2)))
e <- ggplot(data, aes(x = Treatment, y = Netphotosynthesisrate))
e8<-e + geom_boxplot(
  aes(color = Accession), width = 0.5, size = 0.4,
  position = position_dodge(0.8)
  )  +geom_dotplot(
    aes(fill = Accession, color = Accession), binwidth = scale/80, trim = FALSE,
    binaxis='y', stackdir='center', dotsize =1.2,
    position = position_dodge(0.8) )+
  scale_fill_manual(values = c("#F8766D","#619CFF"))+
  scale_color_manual(values = c( "#F8766D","#619CFF"))+ylab(expression(plain(paste("amino acid",~'(','nmol /',~'mg DW',')')))) + theme(text = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+ stat_summary(fun=median, geom='line',aes(group=Accession,col= Accession),position=pd)+ylim(yscale1, yscale2)+annotate("text",annotation$x,annotation$y,label= annotation$label1,parse=TRUE,color="black",size=6/.pt,vjust=0)+annotate("text",annotation$x,annotation$y,label=annotation$label2,parse=TRUE,color="black",size=6/.pt,vjust=1.5)+annotate("text",annotation$x,annotation$y,label= annotation$label3,parse=TRUE,color="black",size=6/.pt,vjust=3)

pdf(paste0("~/Downloads/supplemental_figure2_sheath.pdf"), width=16, height = 9)

ggarrange(e1,e2,e3,e4,e5,e6,e7,e8,labels=c("a","b","c","d","e","f","g","h"),ncol=4,nrow=2, common.legend=TRUE, legend="bottom")
dev.off()

