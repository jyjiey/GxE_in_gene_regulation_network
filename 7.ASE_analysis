
###########
setwd("/Users/jiey/Downloads/")
library(DESeq2)
sample_real_id_matrix<-read.table("/Users/jiey/Downloads/ase.data.csv",header=T,sep=",")
A<- read.table("06_combined_quantified.gene.counts",header=FALSE,sep="\t")
count_matrix<- matrix(0,dim(A)[1],48)
rownames(count_matrix)<-A$V1
colnames(count_matrix)<-c(sample_real_id_matrix$V2)

file.names <- dir("/Users/jiey/Downloads/",pattern ="_combined_quantified.gene.counts")
name<-matrix(0,1,24)

for (i in 1:length(file.names)) {	
matrix <- read.table(file.names[i],header= FALSE,sep="\t")
count_matrix[,i]<- matrix$V2
count_matrix[,c(24+i)]<- matrix$V3
name[,i]<-file.names[i]
}
###########
length(unique(substr(rownames(count_matrix),1,12)))
write.csv(count_matrix,"~/Downloads/ASE_count_matrix.csv")
count_matrix <-read.csv("~/Downloads/ASE_count_matrix.csv", header=T,row.names=1)

Diff<-unique(rownames(which(abs((count_matrix[,1:24]-count_matrix[,25:48]))>0,arr.ind=TRUE)))
Diff_old<-Diff

#newDiff1 <-rownames(count_matrix) [which((rowSums(ifelse((count_matrix[,c(c(6,8,16,2,15,18)+24)]+0.000000001)/(count_matrix[,c(6,8,16,2,15,18)]+ 0.0000001)<0.05,1,0),na.rm=TRUE)>5)& (rowSums(ifelse((count_matrix[,c(7,9,20,11,13,22)]+ 0.0000001)/(count_matrix[,c(c(7,9,20,11,13,22)+24)]+ 0.000000001)<0.05,1,0),na.rm=TRUE)>5))]
#newDiff1 <-rownames(count_matrix)[which((rowSums(ifelse((count_matrix[,c(c(6,8,16,2,15,18)+24)]+0.05)/(count_matrix[,c(6,8,16,2,15,18)]+ 0.05)<0.05,1,0),na.rm=TRUE)>5)& (rowSums(ifelse((count_matrix[,c(7,9,20,11,13,22)]+ 0.05)/(count_matrix[,c(c(7,9,20,11,13,22)+24)]+ 0.05)<0.05,1,0),na.rm=TRUE)>5))]
length(newDiff1)
newDiff1 <-rownames(count_matrix)[which((rowSums(ifelse((count_matrix[,c(c(6,8,16,2,18)+24)]+0.05)/(count_matrix[,c(6,8,16,2,18)]+ 0.05)<0.2,1,0),na.rm=TRUE)>4)& (rowSums(ifelse((count_matrix[,c(7,9,11,13,22)]+ 0.05)/(count_matrix[,c(c(7,9,11,13,22)+24)]+ 0.05)<0.2,1,0),na.rm=TRUE)>4))]


#newDiff1 <-rownames(count_matrix)[which((rowSums(ifelse((count_matrix[,c(c(6,8,16,2,15,18)+24)]+0.05)/(count_matrix[,c(6,8,16,2,15,18)]+ 0.05)<0.05,1,0),na.rm=TRUE)>5)& (rowSums(ifelse((count_matrix[,c(7,9,20,11,13,22)]+ 0.05)/(count_matrix[,c(c(7,9,20,11,13,22)+24)]+ 0.05)<0.05,1,0),na.rm=TRUE)>5))][c(1:1000)]

#ifelse((count_matrix[,c(c(6,8,16,2,15,18)+24)]+0.0000001)/(count_matrix[,c(6,8,16,2,15,18)]+ 0.0000001)<0.05,1,0)

#ifelse((count_matrix[which(rownames(count_matrix)=="Bradi1g06870.4"),c(7,9,20,11,13,22)]+ 0.0000001)/(count_matrix[which(rownames(count_matrix)=="Bradi1g06870.4"),c(c(7,9,20,11,13,22)+24)]+ 0.0000001)<0.05,1,0)count_matrix[which(rownames(count_matrix)=="Bradi1g06870.4"),]
####<0.05, >4 samples
####consider 0?
Diff <-unique(newDiff1[which(newDiff1%in%Diff)])
#87 too few reads
#count_matrix[which(rownames(count_matrix)=="Bradi3g00890"),]
############################
##count_matrix_combined<-count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE")[which(rowSums(count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)]+count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)+24]) >= 5)],c(1:19,21:24)]+count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE")[which(rowSums(count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)]+count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)+24]) >= 5)],c(1:19,21:24)+24]
count_matrix_combined<-count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE")[which(rowSums(count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)]+count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)+24]) > 10)],c(1:19,21:24)]+count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE")[which(rowSums(count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)]+count_matrix[which(rownames(count_matrix)%in% Diff =="FALSE"),c(1:19,21:24)+24]) >10)],c(1:19,21:24)+24]
A<-dim(count_matrix_combined)[1]
count_matrix_combined<-rbind(count_matrix_combined,count_matrix[rownames(count_matrix)%in% Diff,c(1:19,21:24)])
B<-dim(count_matrix_combined)[1]
B1<-count_matrix[rownames(count_matrix)%in% Diff,c(1:19,21:24)+24]
colnames(B1)<-colnames(count_matrix_combined)
count_matrix_combined<-rbind(count_matrix_combined,B1)

C<-dim(count_matrix_combined)[1]

library(DESeq2)
count_matrix_combined_round<-round(count_matrix_combined)

colnames(sample_real_id_matrix)[5]<-"Verables"

#sample_real_id_matrix1<-sample_real_id_matrix[1:24,1:5]
sample_real_id_matrix1<-sample_real_id_matrix[c(1:19,21:24),1:5]
design = formula(~ Verables)
dds <- DESeqDataSetFromMatrix(countData = count_matrix_combined_round,
                              colData = sample_real_id_matrix1,
                              design = design)
dds <- estimateSizeFactors(dds) 
 normalized_counts <- counts(dds, normalized=TRUE)
 
 
 library(factoextra)
colnames(datExpr)
dta<-prcomp(t(normalized_counts))
dta$x[,1]
 pdf("~/Downloads/plot_1.pdf", width=7, height = 5)
fviz_eig(dta, addlabels = TRUE, ylim = c(0, 50))
type<-sample_real_id_matrix1$Verables
type[which(type=="bd3-1")]<-"bd31c"
type[which(type=="droug")]<-"hybrid_d"
type[which(type=="contr")]<-"hybrid_c"
fviz_pca_ind(dta,
             col.ind =type, # Color by the quality of representation
 geom = c("point"),
             repel = TRUE     # Avoid text overlapping
             )
plot(dta$x[,3],dta$x[,4],col= as.factor(type))        
plot(dta$x[,5],dta$x[,6],col= as.factor(type))
plot(dta$x[,7],dta$x[,8],col= as.factor(type))
plot(dta$x[,9],dta$x[,10],col= as.factor(type))     
dev.off()

##bias check


plot(log(rowMeans(normalized_counts[c(A+1):B,])),log(rowMeans(normalized_counts[c(B+1):C,])), pch='.',xlab="Bd21 ASE Counts (Log10)",ylab="Bd3-1 ASE Counts (Log10)",xlim=c(0,12),ylim=c(0,12))
abline(a=0, b=1, col="red")

dim(normalized_counts)
#Count<-rbind(normalized_counts[1:A,c(6,8,16,2,15,18,7,9,11,13)],cbind(normalized_counts[c(A+1):B,c(6,8,16,2,15,18)],normalized_counts[(B+1):C,c(7,9,11,13)]))
#normalized_counts_orignal_unnormalized<-round(cbind(count_matrix_combined_round[c(A+1):B,], count_matrix_combined_round[c(B+1):C,]))

normalized_counts_orignal<-round(cbind(normalized_counts[c(A+1):B,],normalized_counts[c(B+1):C,]))
matrixbd<-read.csv("~/Downloads/matrixbd_for_groupwiseWGCNA.csv", header=T,row.names=1)
#normalized_counts_orignal <-normalized_counts_orignal[which((substr(rownames(normalized_counts_orignal),1,12)%in%substr(colnames(matrixbd),1,12))),]
#sample_real_id_matrix2<-sample_real_id_matrix[-which(sample_real_id_matrix $V2%in%c(87,82,73)),] 
#sample_real_id_matrix2 <-sample_real_id_matrix2[which((sample_real_id_matrix2 $Keep=="1")),]
#sample_real_id_matrix1<-sample_real_id_matrix[-which(sample_real_id_matrix $V2%in%c(87)),]
#normalized_counts_orignal <-normalized_counts_orignal[,-which((sample_real_id_matrix1 $V2%in%c(82,73)))]
#sample_real_id_matrix1<-sample_real_id_matrix1[-which(sample_real_id_matrix1 $V2%in%c(82,73)),]
#normalized_counts_orignal <-normalized_counts_orignal[,which((sample_real_id_matrix1 $Keep==1))]



##remove 82,73
sample_real_id_matrix2<-sample_real_id_matrix[-which(sample_real_id_matrix $V2%in%c(87,82,73)),] 
sample_real_id_matrix2 <-sample_real_id_matrix2[which((sample_real_id_matrix2 $Keep=="1")),]
sample_real_id_matrix1<-sample_real_id_matrix[-which(sample_real_id_matrix $V2%in%c(87)),]
normalized_counts_orignal <-normalized_counts_orignal[,-which((sample_real_id_matrix1 $V2%in%c(82,73)))]
sample_real_id_matrix1<-sample_real_id_matrix1[-which(sample_real_id_matrix1 $V2%in%c(82,73)),]
normalized_counts_orignal <-normalized_counts_orignal[,which((sample_real_id_matrix1 $Keep==1))]




##cis control: hybrid difference; all sample difference
##cis drought: hybrid difference; all sample difference
##cis all: hybrid only; all sample
##cis by E:  hybrid only; all sample 


#normalized_counts_orignal[,which(sample_real_id_matrix2[,8]==1)]<-round((normalized_counts_orignal/2)[,which(sample_real_id_matrix2[,8]==1)])
##hybrid_control
sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,7]=="control")&(sample_real_id_matrix2[,8]=="0")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,7]=="control")&(sample_real_id_matrix2[,8]=="0"))]
design = formula(~  Accession)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
Atable<-as.matrix(results(dds, name='Accession_bd31_vs_bd21'))
cis_c<-ifelse(Atable[,5]<0.1,1,0)
cis_c_q<-ifelse(Atable[,6]<0.1,1,0)
cis_c[which(is.na(cis_c)==TRUE)]<-0
cis_c_q[which(is.na(cis_c_q)==TRUE)]<-0

##hybrid_control
sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,7]=="control")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,7]=="control"))]
design = formula(~  Accession*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
Atable<-as.matrix(results(dds, name='Accession_bd31_vs_bd21'))
cis_c_2<-ifelse(Atable[,5]<0.1,1,0)
cis_c_q_2<-ifelse(Atable[,6]<0.1,1,0)
cis_c_2[which(is.na(cis_c_2)==TRUE)]<-0
cis_c_q_2[which(is.na(cis_c_q_2)==TRUE)]<-0
plyr::count(cbind(cis_c_2, cis_c))


# ##hybrid_control
# Treatment_control<-ifelse(sample_real_id_matrix2$Treatment=="control",0,1)
# sample_real_id_matrix3<-cbind(sample_real_id_matrix2, Treatment_control)
# design = formula(~  Accession*generation*Treatment_control)
# dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              # colData = sample_real_id_matrix3,
                              # design = design)
# "sizeFactors"(dds) <- 1
# dds <- DESeq(dds)
# resultsNames(dds)
# Atable<-as.matrix(results(dds, name='Accession_bd31_vs_bd21'))
# cis_c_3<-ifelse(Atable[,5]<0.1,1,0)
# cis_c_q_3<-ifelse(Atable[,6]<0.1,1,0)
# cis_c_3[which(is.na(cis_c_3)==TRUE)]<-0
# cis_c_q_3[which(is.na(cis_c_q_3)==TRUE)]<-0
# plyr::count(cbind(cis_c_3 ,cis_c_2, cis_c))


sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,7]=="drought")&(sample_real_id_matrix2[,8]=="0")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,7]=="drought")&(sample_real_id_matrix2[,8]=="0"))]
design = formula(~  Accession)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
Atable<-as.matrix(results(dds, name='Accession_bd31_vs_bd21'))
cis_d<-ifelse(Atable[,5]<0.1,1,0)
cis_d_q<-ifelse(Atable[,6]<0.1,1,0)
cis_d[which(is.na(cis_d)==TRUE)]<-0
cis_d_q[which(is.na(cis_d_q)==TRUE)]<-0


##hybrid_control
sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,7]=="drought")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,7]=="drought"))]
design = formula(~  Accession*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
Atable<-as.matrix(results(dds, name='Accession_bd31_vs_bd21'))
cis_d_2<-ifelse(Atable[,5]<0.1,1,0)
cis_d_q_2<-ifelse(Atable[,6]<0.1,1,0)
cis_d_2[which(is.na(cis_d_2)==TRUE)]<-0
cis_d_q_2[which(is.na(cis_d_q_2)==TRUE)]<-0
plyr::count(cbind(cis_d_2, cis_d))


# # ##hybrid_control
# Treatment_control<-ifelse(sample_real_id_matrix2$Treatment=="drought",0,1)
# sample_real_id_matrix3<-cbind(sample_real_id_matrix2, Treatment_control)
# design = formula(~  Accession*generation*Treatment_control)
# dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              # colData = sample_real_id_matrix3,
                              # design = design)
# "sizeFactors"(dds) <- 1
# dds <- DESeq(dds)
# resultsNames(dds)
# Atable<-as.matrix(results(dds, name='Accession_bd31_vs_bd21'))
# cis_d_3<-ifelse(Atable[,5]<0.1,1,0)
# cis_d_q_3<-ifelse(Atable[,6]<0.1,1,0)
# cis_d_3[which(is.na(cis_d_3)==TRUE)]<-0
# cis_d_q_3[which(is.na(cis_d_q_3)==TRUE)]<-0
# plyr::count(cbind(cis_d_3 ,cis_d_2, cis_d))


sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,8]=="0")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,8]=="0"))]
design = formula(~  Accession*Treatment)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
cisxe<-ifelse(as.matrix(results(dds, name="Accessionbd31.Treatmentdrought"))[,5]<0.1,1,0)
cisxe_q<-ifelse(as.matrix(results(dds, name="Accessionbd31.Treatmentdrought"))[,6]<0.1,1,0)
cisxe[which(is.na(cisxe)==TRUE)]<-0
cisxe_q[which(is.na(cisxe_q)==TRUE)]<-0


design = formula(~  Accession*Treatment*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              colData = sample_real_id_matrix2,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
cisxe_2<-ifelse(as.matrix(results(dds, name="Accessionbd31.Treatmentdrought"))[,5]<0.1,1,0)
cisxe_q_2<-ifelse(as.matrix(results(dds, name="Accessionbd31.Treatmentdrought"))[,6]<0.1,1,0)
cisxe_2[which(is.na(cisxe_2)==TRUE)]<-0
cisxe_q_2[which(is.na(cisxe_q_2)==TRUE)]<-0
plyr::count(cbind(cisxe_2, cisxe))



sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,8]=="0")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,8]=="0"))]
design = formula(~  Accession)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
cis<-ifelse(as.matrix(results(dds, name="Accession_bd31_vs_bd21"))[,5]<0.1,1,0)
cis_q<-ifelse(as.matrix(results(dds, name="Accession_bd31_vs_bd21"))[,6]<0.1,1,0)
cis[which(is.na(cis)==TRUE)]<-0
cis_q[which(is.na(cis_q)==TRUE)]<-0


design = formula(~  Accession*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              colData = sample_real_id_matrix2,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
cis_2 <-ifelse(as.matrix(results(dds, name="Accession_bd31_vs_bd21"))[,5]<0.1,1,0)
cis_q_2<-ifelse(as.matrix(results(dds, name="Accession_bd31_vs_bd21"))[,6]<0.1,1,0)
cis_2[which(is.na(cis_2)==TRUE)]<-0
cis_q_2[which(is.na(cis_q_2)==TRUE)]<-0
plyr::count(cbind(cis_2, cis))




cis<-ifelse((cis_2+ cis)>0,1,0)
cis_q<-ifelse((cis_q_2+ cis_q)>0,1,0)
cisxe<-ifelse((cisxe_2+ cisxe)>0,1,0)
cisxe_q<-ifelse((cisxe_q_2+ cisxe_q)>0,1,0)
# cis_c<-ifelse((cis_c_3+cis_c_2+ cis_c)>0,1,0)
# cis_c_q<-ifelse((cis_c_q_3+cis_c_q_2+ cis_c_q)>0,1,0)
# cis_d<-ifelse((cis_d_3+cis_d_2+ cis_d)>0,1,0)
# cis_d_q<-ifelse((cis_d_q_3+cis_d_q_2+ cis_d_q)>0,1,0)
cis_c<-ifelse((cis_c_2+ cis_c)>0,1,0)
cis_c_q<-ifelse((cis_c_q_2+ cis_c_q)>0,1,0)
cis_d<-ifelse((cis_d_2+ cis_d)>0,1,0)
cis_d_q<-ifelse((cis_d_q_2+ cis_d_q)>0,1,0)

cis_both<-ifelse((cis_c+ cis_d)>1,1,0)
cis_both_q<-ifelse((cis_c_q+ cis_d_q)>1,1,0)
cis_e<-ifelse((cis_c+ cis_d)==1,1,0)
cis_e_q<-ifelse((cis_c_q+ cis_d_q)==1,1,0)

##trans c 
##trans d
##trans all
##trans by E 

##hybrid_control
sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,7]=="drought")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,7]=="drought"))]
design = formula(~  Accession*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
Atable<-as.matrix(results(dds, name='Accessionbd31.generation'))
trans_d<-ifelse(Atable[,5]<0.1,1,0)
trans_d_q<-ifelse(Atable[,6]<0.1,1,0)
trans_d[which(is.na(trans_d)==TRUE)]<-0
trans_d_q[which(is.na(trans_d_q)==TRUE)]<-0


# ##hybrid_control
# Treatment_control<-ifelse(sample_real_id_matrix2$Treatment=="drought",0,1)
# sample_real_id_matrix3<-cbind(sample_real_id_matrix2, Treatment_control)
# design = formula(~  Accession*generation* Treatment_control)
# dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              # colData = sample_real_id_matrix3,
                              # design = design)
# "sizeFactors"(dds) <- 1
# dds <- DESeq(dds)
# resultsNames(dds)
# Atable<-as.matrix(results(dds, name='Accessionbd31.generation'))
# trans_d_2<-ifelse(Atable[,5]<0.1,1,0)
# trans_d_q_2<-ifelse(Atable[,6]<0.1,1,0)
# trans_d_2[which(is.na(trans_d_2)==TRUE)]<-0
# trans_d_q_2[which(is.na(trans_d_q_2)==TRUE)]<-0




##hybrid_control
sample_real_id_matrix3 <-sample_real_id_matrix2[which((sample_real_id_matrix2[,7]=="control")),]
normalized_counts_orignal_parent_control<-normalized_counts_orignal[,which((sample_real_id_matrix2[,7]=="control"))]
design = formula(~  Accession*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal_parent_control,
                              colData = sample_real_id_matrix3,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
resultsNames(dds)
Atable<-as.matrix(results(dds, name='Accessionbd31.generation'))
trans_c<-ifelse(Atable[,5]<0.1,1,0)
trans_c_q<-ifelse(Atable[,6]<0.1,1,0)
trans_c[which(is.na(trans_c)==TRUE)]<-0
trans_c_q[which(is.na(trans_c_q)==TRUE)]<-0

# ##hybrid_control
# Treatment_control<-ifelse(sample_real_id_matrix2$Treatment=="control",0,1)
# sample_real_id_matrix3<-cbind(sample_real_id_matrix2, Treatment_control)
# design = formula(~  Accession*generation* Treatment_control)
# dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              # colData = sample_real_id_matrix3,
                              # design = design)
# "sizeFactors"(dds) <- 1
# dds <- DESeq(dds)
# resultsNames(dds)
# Atable<-as.matrix(results(dds, name='Accessionbd31.generation'))
# trans_c_2<-ifelse(Atable[,5]<0.1,1,0)
# trans_c_q_2<-ifelse(Atable[,6]<0.1,1,0)
# trans_c_2[which(is.na(trans_c_2)==TRUE)]<-0
# trans_c_q_2[which(is.na(trans_c_q_2)==TRUE)]<-0


design = formula(~  Accession*Treatment*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              colData = sample_real_id_matrix2,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
transxe<-ifelse(as.matrix(results(dds, name="Accessionbd31.Treatmentdrought.generation"))[,5]<0.1,1,0)
transxe_q<-ifelse(as.matrix(results(dds, name="Accessionbd31.Treatmentdrought.generation"))[,6]<0.1,1,0)
transxe[which(is.na(transxe)==TRUE)]<-0
transxe_q[which(is.na(transxe_q)==TRUE)]<-0


design = formula(~  Accession*generation)
dds <- DESeqDataSetFromMatrix(countData = normalized_counts_orignal,
                              colData = sample_real_id_matrix2,
                              design = design)
"sizeFactors"(dds) <- 1
dds <- DESeq(dds)
trans<-ifelse(as.matrix(results(dds, name="Accessionbd31.generation"))[,5]<0.1,1,0)
trans_q<-ifelse(as.matrix(results(dds, name="Accessionbd31.generation"))[,6]<0.1,1,0)
trans[which(is.na(trans)==TRUE)]<-0
trans_q[which(is.na(trans_q)==TRUE)]<-0

# # cbind(trans_c_2, trans_c, trans_c_q_2, trans_c_q, trans_d_2, trans_d, trans_d_q_2, trans_d_q)[which(names(cis_2)%in%c("Bradi3g00920", "Bradi1g51770", "Bradi1g15840", "Bradi1g53680", "Bradi3g45230")),]
# trans_c<-ifelse((trans_c_2+ trans_c)>0,1,0)
# trans_c_q<-ifelse((trans_c_q_2+ trans_c_q)>0,1,0)
# trans_d<-ifelse((trans_d_2+ trans_d)>0,1,0)
# trans_d_q<-ifelse((trans_d_q_2+ trans_d_q)>0,1,0)


trans_c<-ifelse(( trans_c)>0,1,0)
trans_c_q<-ifelse(( trans_c_q)>0,1,0)
trans_d<-ifelse(( trans_d)>0,1,0)
trans_d_q<-ifelse(( trans_d_q)>0,1,0)

trans_both<-ifelse((trans_c + trans_d)>1,1,0)
trans_both_q<-ifelse((trans_c_q+ trans_d_q)>1,1,0)
trans_e<-ifelse((trans_c+ trans_d)==1,1,0)
trans_e_q<-ifelse((trans_c_q+ trans_d_q)==1,1,0)



 diff_ge <-read.csv("~/Desktop/draft_SI/draft_data/gxe_extended.csv", head=TRUE)
 ge <-diff_ge[,1]
 diff_g <-read.csv("~/Desktop/draft_SI/draft_data/g_extended.csv", head=TRUE)
 g<-diff_g[,1]
 diff_g <-read.csv("~/Desktop/draft_SI/draft_data/e_extended.csv", head=TRUE)
 e<-diff_g[,1]
dge_gene<-rep(0,dim(Atable)[1])
dge_gene[which(substr(rownames(Atable),1,12)%in%substr(g,1,12))]<-"g"
dge_gene[which(substr(rownames(Atable),1,12)%in%substr(e,1,12))]<-"e"
dge_gene[which((substr(rownames(Atable),1,12)%in%substr(g,1,12))&(substr(rownames(Atable),1,12)%in%substr(e,1,12)))]<-"g+e"
dge_gene[which((substr(rownames(Atable),1,12)%in%substr(ge,1,12)))]<-"gxe"

data<-cbind(cis,cisxe,cis_both,cis_e,trans,transxe,trans_both,trans_e)
colnames(data)<-c("cis","cisxe","indirect_cis_both","indirect_cis_e","trans","transxe","indirect_trans_both","indirect_trans_e")
data_ase<-data
write.csv(data_ase,"~/Downloads/data_ase.csv")

none<-ifelse(rowSums(data)==0,1,0)
data<-cbind(data,none)
data_q<-cbind(cis_q,cisxe_q,cis_both_q,cis_e_q,trans_q,transxe_q,trans_both_q,trans_e_q)
colnames(data_q)<-c("cis","cisxe","indirect_cis_both","indirect_cis_e","trans","transxe","indirect_trans_both","indirect_trans_e")
none_q<-ifelse(rowSums(data_q)==0,1,0)
data_q<-cbind(data_q,none_q)
data_q<-data_q[,-which(colSums(data_q)==0)]

data_ase<-data


data_ase [which(rownames(data_ase)%in%c("Bradi3g00920", "Bradi1g51770", "Bradi1g15840", "Bradi1g53680", "Bradi3g45230")),]
#library(ggplot2)
 pdf("~/Downloads/plot3.pdf", width=3, height = 3.1)
 a<-c("Bradi3g00920", "Bradi1g51770", "Bradi1g15840", "Bradi1g53680", "Bradi3g45230")
for(i in 1:length(a)){
	names1<-Names<-substr(a,1,12)[i]
	#names1<-Names<-"Bradi3g45230"
data<-cbind(sample_real_id_matrix2$Accession, sample_real_id_matrix2 $Treatment, sample_real_id_matrix2 $generation, c(normalized_counts_orignal[which(rownames(normalized_counts_orignal)== Names),]))
colnames(data)<-c("Accession","Treatment","generation","expression")
write.csv(data,"~/Data_temperatl.csv")
data <-read.csv("~/Data_temperatl.csv", header=T)
boxplot(expression ~ Treatment*Accession*generation,data=data,xlab=paste0(cis[which(dge_gene=="gxe")[i]],cisxe[which(dge_gene=="gxe")[i]],cis_both[which(dge_gene=="gxe")[i]],cis_e[which(dge_gene=="gxe")[i]],trans[which(dge_gene=="gxe")[i]],transxe[which(dge_gene=="gxe")[i]],trans_both[which(dge_gene=="gxe")[i]],trans_e[which(dge_gene=="gxe")[i]]))

if(length(which(substr(colnames(matrixbd),1,12)== Names))>0){

plot(matrixbd[,which(substr(colnames(matrixbd),1,12)== Names)])

groupg<-c(rep(0,48),rep(0,48),rep(1,48),rep(1,48))
groupe<-c(rep(0,48),rep(1,48),rep(0,48),rep(1,48))
boxplot(matrixbd[,which(substr(colnames(matrixbd),1,12)== Names)]~groupe*groupg)
}


ggplot(data, aes(Treatment, expression,group= interaction(Treatment, Accession,generation),colour=interaction(Treatment,genotype) ))  + geom_boxplot(fill = c(rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1),rep("blue",1),rep("red", 1),rep("purple", 1),rep("orange", 1)), color = c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(
    fun.y = median,
    geom = 'line',
   aes(Treatment, expression,group=interaction(Accession, generation)),color=c(rep("black",1),rep("black", 1),rep("black", 1),rep("black", 1),rep("grey",1),rep("grey", 1),rep("grey", 1),rep("grey", 1)),
    position = position_dodge(width = 0.8) #this has to be added
  )+ylim(min(data$expression)-0.5,max(data$expression)+2)+annotate("text",1,max(data$expression),label="F1",color="black",size=5,vjust=0)+annotate("text",1,max(data$expression)*0.9,label= "F0",color="grey",size= 5,vjust=1.5)+ylab(paste0(Names,"gene expression"))+ theme(
panel.background = element_blank())
}
dev.off()


