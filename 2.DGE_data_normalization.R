##section_3.data for modules
#### normalize and PCA
#RNA_Brachy <-read.csv("/scratch/users/jyun/chapter1/RNAseq/expected_count.csv", header=T,row.names=1)
RNA_Brachy <-read.csv("~/Desktop/all things/expected_count_with_cloroplast_genes.csv", header=T,row.names=1)
setwd("/Users/jiey/Desktop/final")
#setwd("/scratch/users/jyun/chapter1/RNAseq")
length(which(rowMax(RNA_Brachy)>0))
#RNA_Brachy_other <-read.csv("/scratch/users/jyun/chapter1/RNAseq/Data_otherinfo_updatedbatch.csv", header=T)
RNA_Brachy_other <-read.csv("~/Desktop/codes jwafs/Data_otherinfo_updatedbatch.csv", header=T)
head(RNA_Brachy_other)
RNA_Brachy_other<-RNA_Brachy_other[-210,c(2,3,4,16,17,18,19,21,23,24,25,27,28,29,30)]
dim(RNA_Brachy)
RNA_Brachy<-as.matrix(RNA_Brachy[,-210])
library(DESeq2)
RNA_Brachy<-round(RNA_Brachy)
RNA_Brachy_other$Accessiontreatment<-paste0(RNA_Brachy_other$Accession, RNA_Brachy_other$Treatment)
RNA_Brachy_other$Accessiontreatment<-as.factor(RNA_Brachy_other$Accessiontreatment)
dim(RNA_Brachy_other)
#RNA_Brachy<-RNA_Brachy[,-which(colnames(RNA_Brachy)=="X80")]
#RNA_Brachy_other<-RNA_Brachy_other[-which(RNA_Brachy_other[,1]==80),]
design = formula(~ Accession + Treatment + Accession: Treatment)
dds <- DESeqDataSetFromMatrix(countData = RNA_Brachy,
                              colData = RNA_Brachy_other,
                              design = design)                           
 dds <- estimateSizeFactors(dds) 
 sizeFactors(dds)
 normalized_counts <- counts(dds, normalized=TRUE)
 datExpr = log(normalized_counts +1,2)
library(factoextra)
colnames(datExpr)
dta<-prcomp(t(datExpr))
dta$x[,1]
fviz_eig(dta, addlabels = TRUE, ylim = c(0, 50))
type<-rep(0,213)
type[which(RNA_Brachy_other$Accessiontreatment=="BD21drought")]<-"Bd21drought"
type[which(RNA_Brachy_other$Accessiontreatment=="BD21normal")]<-"Bd21control"
type[which(RNA_Brachy_other$Accessiontreatment=="BD3-1drought")]<-"Bd3-1drought"
type[which(RNA_Brachy_other$Accessiontreatment=="BD3-1normal")]<-"Bd3-1control"
fviz_pca_ind(dta,
             col.ind =type, # Color by the quality of representation
 geom = c("point"),
             repel = TRUE     # Avoid text overlapping
             )
b<-fviz_pca_ind(dta,
            col.ind =type, 
 geom = c("point"),
             repel = TRUE    
             )+ scale_color_manual(values=c("red","blue","orange","purple"))




##DGE with >10 samples genes and plot
design = formula(~ Accession + Treatment + Accession: Treatment)
dds <- DESeqDataSetFromMatrix(countData = RNA_Brachy,
                              colData = RNA_Brachy_other,
                              design = design)                              
 dds <- estimateSizeFactors(dds) 
 sizeFactors(dds)
 normalized_counts <- counts(dds, normalized=TRUE)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)
res_gxe<-results(dds,alpha=0.1)
dds1<-dds
res_g <-results(dds,contrast=c('Accession','BD21','BD3-1') , alpha=0.05)
res_e<-results(dds,contrast=c('Treatment','drought','normal') , alpha=0.05)
length(which((res_e[,2]>0)&(res_e[,6]<0.05)))
write.csv(cbind(rownames(res_gxe), res_gxe[,c(2,6)], res_g[,c(2,6)], res_e[,c(2,6)]), "~/Downloads/gxe_lfc.csv")
design2 = formula(~ Accessiontreatment)
dds2 <- DESeqDataSetFromMatrix(countData = RNA_Brachy,
                              colData = RNA_Brachy_other,
                              design = design2)
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]
dds2 <- DESeq(dds2)
res_bd21droughtvsnormal<-results(dds2,contrast=c('Accessiontreatment','BD21drought','BD21normal') , alpha=0.05)
res_bd31droughtvsnormal<-results(dds2,contrast=c('Accessiontreatment','BD3-1drought','BD3-1normal') , alpha=0.05)
#both higher
length(which((res_bd31droughtvsnormal[,2]>0)&(res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,2]>0)&(res_bd21droughtvsnormal[,6]<0.05)&(res_gxe[,6]<0.1)))
#both lower
length(which((res_bd31droughtvsnormal[,2]<0)&(res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,2]<0)&(res_bd21droughtvsnormal[,6]<0.05)&(res_gxe[,6]<0.1)))
#bd31 lower
length(which((res_bd31droughtvsnormal[,2]<0)&(res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,6]>0.05)&(res_gxe[,6]<0.1)))
#bd31 higher
length(which((res_bd31droughtvsnormal[,2]>0)&(res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,6]>0.05)&(res_gxe[,6]<0.1)))
#bd21 higher
length(which((res_bd21droughtvsnormal[,2]>0)&(res_bd21droughtvsnormal[,6]<0.05)&(res_bd31droughtvsnormal[,6]>0.05)&(res_gxe[,6]<0.1)))
#bd21 lower
length(which((res_bd21droughtvsnormal[,2]<0)&(res_bd21droughtvsnormal[,6]<0.05)&(res_bd31droughtvsnormal[,6]>0.05)&(res_gxe[,6]<0.1)))
#opposite
length(which((res_bd31droughtvsnormal[,2]<0)&(res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,2]>0)&(res_bd21droughtvsnormal[,6]<0.05)&(res_gxe[,6]<0.1)))
length(which((res_bd31droughtvsnormal[,2]>0)&(res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,2]<0)&(res_bd21droughtvsnormal[,6]<0.05)&(res_gxe[,6]<0.1)))
#plot Bd21 vs Bd3-1
GxE <-ifelse(res_gxe[,6]<0.1,1,0)
Category<-rep("Bd21 ns; Bd3-1 ns",length(GxE))
Category[which((res_bd31droughtvsnormal[,6]>0.05)&(res_bd21droughtvsnormal[,6]<0.05))]<-"Bd21 E; Bd3-1 ns"
Category[which((res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,6]<0.05))]<-"Bd21 E; Bd3-1 E"
Category[which((res_bd31droughtvsnormal[,6]<0.05)&(res_bd21droughtvsnormal[,6]> 0.05))]<-"Bd21 ns; Bd3-1 E"
bd21<-res_bd21droughtvsnormal[,2]
bd31<-res_bd31droughtvsnormal[,2]
Category[which(GxE==1)]<-"GxE"
matrix<-cbind(bd21, bd31,Category)
cols=c("blue","yellow","green","grey","red")
Data<-as.data.frame(matrix)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)

df2<-matrix[which(rownames(res_bd21droughtvsnormal)%in%c("Bradi3g00920.1","Bradi1g65130.1")),]

write.csv(df2,"~/df2_temperatl.csv")
df2 <-read.csv("~/df2_temperatl.csv", header=T,row.names=1)
#pdf("~/Downloads/plot_dge.pdf", width=7, height = 5)
a<-ggplot(Data, aes(x = bd21, y = bd31, color = Category)) +
  geom_point() +xlab("Lfc of E effect in Bd21")+ylab("Lfc of E effect in Bd3-1")+xlim(-5,5)+ylim(-5,5)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank())+scale_color_manual(values = alpha(cols,0.5))+ geom_hline(yintercept = 0, color = "grey") + geom_vline(xintercept = 0, color = "grey") + geom_abline(yintercept = 0, slope = 1, color = "grey") + 
  geom_point(data= df2, aes(x= bd21, y= bd31), pch=21, 
             fill=NA, size=4, 
             colour="red", stroke=1) + annotate("text", x = c(3.5,-2), y = c(0.5,-5), label = c("Bradi3g00920", "Bradi1g65130"))




###normalization with all genes and remove <0.1 variance genes
RNA_Brachy <-read.csv("~/Desktop/all things/expected_count_with_cloroplast_genes.csv", header=T,row.names=1)
length(which(rowMax(RNA_Brachy)>0))
RNA_Brachy_other <-read.csv("~/Desktop/codes jwafs/Data_otherinfo_updatedbatch.csv", header=T)
head(RNA_Brachy_other)
##sample210 weird and removed 
RNA_Brachy_other<-RNA_Brachy_other[-210,c(2,3,4,16,17,18,19,21,23,24,25,27,28,29,30)]
dim(RNA_Brachy)
RNA_Brachy<-as.matrix(RNA_Brachy[,-210])
library(DESeq2)
RNA_Brachy<-round(RNA_Brachy)
RNA_Brachy_other$Accessiontreatment<-paste0(RNA_Brachy_other$Accession, RNA_Brachy_other$Treatment)
RNA_Brachy_other$Accessiontreatment<-as.factor(RNA_Brachy_other$Accessiontreatment)
##sample80 removed
#RNA_Brachy<-RNA_Brachy[,-which(colnames(RNA_Brachy)=="X80")]
#RNA_Brachy_other<-RNA_Brachy_other[-which(RNA_Brachy_other[,1]==80),]
design = formula(~ Accession + Treatment + Accession: Treatment)
dds <- DESeqDataSetFromMatrix(countData = RNA_Brachy,
                              colData = RNA_Brachy_other,
                              design = design)  
                              keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]                         
dds <- estimateSizeFactors(dds) 
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
datExpr = log(normalized_counts +1,2)

#write.csv(datExpr, "~/Desktop/all things/expected_count_with_cloroplast_genes_normalized.csv")
set.seed(321)
dim(Bd31normal)
Bd31normal<-datExpr[,sample(which((RNA_Brachy_other[,2]=="BD3-1")&(RNA_Brachy_other[,3]=="normal")),48)]
Bd21normal<-datExpr[,sample(which((RNA_Brachy_other[,2]=="BD21")&(RNA_Brachy_other[,3]=="normal")),48)]
Bd31drought<-datExpr[,sample(which((RNA_Brachy_other[,2]=="BD3-1")&(RNA_Brachy_other[,3]=="drought")),48)]
Bd21drought<-datExpr[,sample(which((RNA_Brachy_other[,2]=="BD21")&(RNA_Brachy_other[,3]=="drought")),48)]
matrixbdorigianl<-cbind(Bd21normal, Bd21drought, Bd31normal, Bd31drought)
matrixbd<-cbind(Bd21normal, Bd21drought, Bd31normal, Bd31drought)[-unique(c(which(rowVars(Bd21normal)<0.1),which(rowVars(Bd31normal)<0.1),which(rowVars(Bd31drought)<0.1),which(rowVars(Bd21drought)<0.1))),]
#write.csv(matrixbd, "~/Downloads/matrixbd_groupwise_module.csv")




###gene annotation
l<-rownames(res_gxe)[which(res_gxe[,6]<0.1)]
h<-rep(1,length(l))
if(length(l)>0){
RNA_Brachy_ANNOTATION3<-read.csv("~/Desktop/gene_anno_phytozome.csv", header=T)
for (j in 1: length(l)){
if(l[j]%in% RNA_Brachy_ANNOTATION3[,1]){
h[j]<-as.character(RNA_Brachy_ANNOTATION3[which(RNA_Brachy_ANNOTATION3[,1]==l[j])[1],2])
h[j]<-ifelse("lcl.LT558597"==substr(l[j],1,12),"chloroplast gene",h[j])
}
}
}
plyr::count(h)
h<-rep(1,length(l))
if(length(l)>0){
RNA_Brachy_ANNOTATION3<-read.csv("~/Desktop/codes jwafs/annotation_info_PART3.csv", header=T)
#~/Desktop/codes jwafs/annotation_info_PART3.csv
for (j in 1: length(l)){
if(substr(as.character(l[j]),1,12)%in% RNA_Brachy_ANNOTATION3[,1]){
h[j]<-as.character(RNA_Brachy_ANNOTATION3[which(RNA_Brachy_ANNOTATION3[,1]==substr(l[j],1,12))[1],3])
h[j]<-ifelse("lcl.LT558597"==substr(l[j],1,12),"chloroplast gene",h[j])
}
}
}

#GO enrichment
a<-substr(l,6,nchar(l)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
A<- paste0(BRADI,a,v3)
a<-as.data.frame(A)
write.csv(a, "~/Downloads/a.csv")
#use pather webpage for GO enrichment check



### individual plot 
matrixbdorigianl <-t(matrixbdorigianl)
size1=3
names1<-"Bradi1g65130.1"
value<-as.numeric(matrixbdorigianl[,which(colnames(matrixbdorigianl)==names1)])
type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
genotype<-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
condition<-c(rep("control",48), rep("drought",48), rep("control",48), rep("drought",48))
gene<-rep(1,192)
#pdf(paste0("~/Downloads/plot_one_gene.pdf"), width=3, height = 3)
Data2 <-cbind(condition, genotype, value, type,gene)
Data2<-as.data.frame(Data2)
write.csv(Data2,"~/Data_temperatl.csv")
Data2<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
library(ggplot2)
  d<-ggplot(Data2, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1)),position = position_dodge(width = 0.8), outlier.shape = NA,alpha=0.5) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+annotate("text",2,max(value)-0.5,label=substr(names1,1,12),color="black",size=size1,vjust=0)+ylab("transcript abundance")


#dev.off()
names1<-"Bradi3g00920.1"
value<-as.numeric(matrixbdorigianl[,which(colnames(matrixbdorigianl)==names1)])
type<-c(rep("Bd21c",48), rep("Bd21d",48), rep("Bd3-1c",48), rep("Bd3-1d",48))
genotype<-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
condition<-c(rep("control",48), rep("drought",48), rep("control",48), rep("drought",48))
gene<-rep(1,192)
#pdf(paste0("~/Downloads/plot_one_gene.pdf"), width=3, height = 3)
Data1 <-cbind(condition, genotype, value, type,gene)
Data1 <-as.data.frame(Data1)
write.csv(Data1,"~/Data_temperatl.csv")
Data1<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
library(ggplot2)
c<-ggplot(Data1, aes(condition, value,group= interaction(condition,genotype,gene),colour=interaction(condition,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1)),position = position_dodge(width = 0.8), outlier.shape = NA,alpha=0.5) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(condition, value,group=interaction(genotype,gene)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+annotate("text",1,max(value)-0.5,label=substr(names1,1,12),color="black",size=size1,vjust=0)+ylab("transcript abundance")


pdf("~/Downloads/plot_figure3.pdf", width=10, height =7)
e<-NULL
ggarrange(ggarrange(b,a,labels=c("A","B"),ncol=2,nrow=1),ggarrange(c,e,d,e,labels=c("C","","D",""),common.legend=TRUE, legend="bottom",ncol=4,nrow=1, heights=0.8, widths = c(1,0.5,1,0.5)),ncol=1,nrow=2,align = "v")

dev.off()             
           
