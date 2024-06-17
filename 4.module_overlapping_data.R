##this file is run with two sh file in 81 parallel jobs, for example
#  #!/bin/bash
#  #!/usr/bin/env bash
#  SUBSAMPLE_IDS=({1..81})    
#  SUBSAMPLE_ID="${SUBSAMPLE_IDS[$1]}"
#  module load  R/2020a 
#  Rscript /home/gridsan/jyun/network/parallel_10_13_supercloud_correctedlabel_202312.R ${SUBSAMPLE_ID} > /home/gridsan/jyun/network/parallel_10_13_supercloud_202312.R_${SUBSAMPLE_ID}_real.out
#  exit 0
#
#  #!/bin/bash
#  #!/usr/bin/env bash
#  #SBATCH --nodes=1               # Number of nodes
#  #SBATCH --ntasks-per-node=1          # MPI processes per node
#  #SBATCH -a 0-80
#  ./run_network.sh $SLURM_ARRAY_TASK_ID
#

#library(qvalue)
library(RColorBrewer)
#library(GO.db)
#library(rbioapi)
#library(igraph)
library(plyr)
library(Rcpp)
library(tibble)
#
install.packages("SuperExactTest",repos="https://cran.r-project.org")
library(SuperExactTest)
 packageVersion("SuperExactTest") 

#install.packages("SuperExactTest",repos="https://cran.r-project.org",lib="/scratch/users/jyun/chapter1/RNAseq/Rresults")
#library(SuperExactTest,lib="/scratch/users/jyun/chapter1/RNAseq/Rresults")
#install.packages("SuperExactTest",repos="https://cran.r-project.org")
 a<-as.character(c(1,12,13,14))
b<-as.character(c(5,6,17))
apply(as.matrix(expand.grid(as.character(a),"-",as.character(b),"-")), 1, paste, collapse=" ")
#scratch/users/jyun/chapter1/RNAseq/Rresults/groupwise_WGCNA/
#home/gridsan/jyun/network/
matrixbd<-read.csv("/home/gridsan/jyun/network/matrixbd_for_groupwiseWGCNA.csv", header=T,row.names=1)
cores<-read.csv("/home/gridsan/jyun/network/TOM1_NEW_membership_2_0.99_30_mid_10_1.csv", header=T,row.names=1)
bd21cmembership0 <-cores[,1]
TOM2_NEW_membership<-read.csv("/home/gridsan/jyun/network/TOM2_NEW_membership_2_0.99_30_mid_10_1.csv", header=T,row.names=1)
bd21dmembership0 <-TOM2_NEW_membership[,1]
TOM3_NEW_membership<-read.csv("/home/gridsan/jyun/network/TOM3_NEW_membership_2_0.99_30_mid_10_1.csv", header=T,row.names=1)
bd31cmembership0 <-TOM3_NEW_membership[,1]
TOM4_NEW_membership<-read.csv("/home/gridsan/jyun/network/TOM4_old_smaller.csv", header=T,row.names=1)
bd31dmembership0<-TOM4_NEW_membership[,1]
TOM1_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM1_NEW_membership_3_0.984_15_break_10_1.csv", header=T,row.names=1)
bd21cmembershiploose <-TOM1_NEW_membership_para4[,1]
TOM2_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM2_NEW_membership_3_0.984_15_break_10_1.csv", header=T,row.names=1)
bd21dmembershiploose <-TOM2_NEW_membership_para4[,1]
TOM3_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM3_NEW_membership_3_0.984_15_break_10_1.csv", header=T,row.names=1)
bd31cmembershiploose <-TOM3_NEW_membership_para4[,1]
TOM4_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM4_old_smaller2.csv", header=T,row.names=1)
bd31dmembershiploose <-TOM4_NEW_membership_para4[,1]
j=1
for (i in as.numeric(levels(as.factor(bd21cmembershiploose))[-1])){
bd21cmembershiploose[which(bd21cmembershiploose==i, arr.ind = T)]<-j
	j=j+1
}
j=1
for (i in as.numeric(levels(as.factor(bd21dmembershiploose))[-1])){
bd21dmembershiploose[which(bd21dmembershiploose==i, arr.ind = T)]<-j
	j=j+1
}
j=1
for (i in as.numeric(levels(as.factor(bd31cmembershiploose))[-1])){
bd31cmembershiploose[which(bd31cmembershiploose==i, arr.ind = T)]<-j
	j=j+1
}

j=1
for (i in as.numeric(levels(as.factor(bd31dmembershiploose))[-1])){
bd31dmembershiploose[which(bd31dmembershiploose==i, arr.ind = T)]<-j
	j=j+1
}

TOM1_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM1_NEW_membership_1_0.99_30_big_10_1.csv", header=T,row.names=1)
bd21cmembershipstrict <-TOM1_NEW_membership_para4[,1]
TOM2_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM2_NEW_membership_1_0.99_30_big_10_1.csv", header=T,row.names=1)
bd21dmembershipstrict <-TOM2_NEW_membership_para4[,1]
TOM3_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM3_NEW_membership_1_0.99_30_big_10_1.csv", header=T,row.names=1)
bd31cmembershipstrict <-TOM3_NEW_membership_para4[,1]
TOM4_NEW_membership_para4<-read.csv("/home/gridsan/jyun/network/TOM4_old_smaller3.csv", header=T,row.names=1)
bd31dmembershipstrict <-TOM4_NEW_membership_para4[,1]

bd21cmembership<-cbind(bd21cmembership0,bd21cmembershiploose,bd21cmembershipstrict)
bd21dmembership<-cbind(bd21dmembership0,bd21dmembershiploose,bd21dmembershipstrict)
bd31cmembership<-cbind(bd31cmembership0,bd31cmembershiploose,bd31cmembershipstrict)
bd31dmembership<-cbind(bd31dmembership0,bd31dmembershiploose,bd31dmembershipstrict)

##all change to loose network module number: check is it the same for reference module number to check if there's any difference 
total=16207
pcutoff=10e-12
bd21cmembership<-cbind(bd21cmembership0,bd21cmembershiploose,bd21cmembershipstrict)

TOM1_NEW_membership_para1<-bd21cmembership[,2]
TOM1_NEW_membership_para4<-bd21cmembership[,1]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21cmembership<-cbind(bd21cmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21cmembership<-cbind(bd21cmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21cmembership<-cbind(bd21cmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd21cmembership[,2])
bd21cmembership_mid_breakdown_list<-bd21cmembership_mid_unique <-bd21cmembership_break_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd21cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd21cmembership[which(bd21cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd21cmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}
max_number
###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd21cmembership_mid_breakdown_list<-c(bd21cmembership_mid_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd21cmembership_break_unique <-levels(as.factor(bd21cmembership[,2]))[-which(levels(as.factor(bd21cmembership[,2]))%in%unique(as.vector(replacement)))]
###add mid unique name to 0 

if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd21cmembership[which(bd21cmembership[,1]==levels(as.factor(bd21cmembership[,1]))[-1][which(colSums(overlap)==0)]),4:6]<-0
	}
}


TOM1_NEW_membership_para1<-bd21cmembership[,2]
TOM1_NEW_membership_para4<-bd21cmembership[,3]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21cmembership<-cbind(bd21cmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21cmembership<-cbind(bd21cmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21cmembership<-cbind(bd21cmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd21cmembership[,2])
bd21cmembership_bigger_breakdown_list<-bd21cmembership_bigger_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd21cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd21cmembership[which(bd21cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd21cmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}


###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd21cmembership_bigger_breakdown_list<-c(bd21cmembership_bigger_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd21cmembership_break_unique <-c(bd21cmembership_break_unique ,levels(as.factor(bd21cmembership[,2]))[-which(levels(as.factor(bd21cmembership[,2]))%in%unique(as.vector(replacement)))])
###add bigger unique name to all at the end and bd21cmembership_bigger_unique
if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd21cmembership[which(bd21cmembership[,1]==levels(as.factor(bd21cmembership[,1]))[-1][which(colSums(overlap)==0)]),7:9]<-0
	}
}


############################################
bd21dmembership<-cbind(bd21dmembership0,bd21dmembershiploose,bd21dmembershipstrict)

TOM1_NEW_membership_para1<-bd21dmembership[,2]
TOM1_NEW_membership_para4<-bd21dmembership[,1]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21dmembership<-cbind(bd21dmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21dmembership<-cbind(bd21dmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21dmembership<-cbind(bd21dmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd21dmembership[,2])
bd21dmembership_mid_breakdown_list<-bd21dmembership_mid_unique <-bd21dmembership_break_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd21dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd21dmembership[which(bd21dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd21dmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}
max_number
###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd21dmembership_mid_breakdown_list<-c(bd21dmembership_mid_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd21dmembership_break_unique <-levels(as.factor(bd21dmembership[,2]))[-which(levels(as.factor(bd21dmembership[,2]))%in%unique(as.vector(replacement)))]
###add mid unique name to 0 

if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd21dmembership[which(bd21dmembership[,1]==levels(as.factor(bd21dmembership[,1]))[-1][which(colSums(overlap)==0)]),4:6]<-0
	}
}


TOM1_NEW_membership_para1<-bd21dmembership[,2]
TOM1_NEW_membership_para4<-bd21dmembership[,3]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21dmembership<-cbind(bd21dmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21dmembership<-cbind(bd21dmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd21dmembership<-cbind(bd21dmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd21dmembership[,2])
bd21dmembership_bigger_breakdown_list<-bd21dmembership_bigger_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd21dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd21dmembership[which(bd21dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd21dmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}


###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd21dmembership_bigger_breakdown_list<-c(bd21dmembership_bigger_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd21dmembership_break_unique <-c(bd21dmembership_break_unique ,levels(as.factor(bd21dmembership[,2]))[-which(levels(as.factor(bd21dmembership[,2]))%in%unique(as.vector(replacement)))])
###add bigger unique name to all at the end and bd21dmembership_bigger_unique
if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd21dmembership[which(bd21dmembership[,1]==levels(as.factor(bd21dmembership[,1]))[-1][which(colSums(overlap)==0)]),7:9]<-0
	}
}


########################
bd31cmembership<-cbind(bd31cmembership0,bd31cmembershiploose,bd31cmembershipstrict)

TOM1_NEW_membership_para1<-bd31cmembership[,2]
TOM1_NEW_membership_para4<-bd31cmembership[,1]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31cmembership<-cbind(bd31cmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31cmembership<-cbind(bd31cmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31cmembership<-cbind(bd31cmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd31cmembership[,2])
bd31cmembership_mid_breakdown_list<-bd31cmembership_mid_unique <-bd31cmembership_break_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd31cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd31cmembership[which(bd31cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd31cmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}
max_number
###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd31cmembership_mid_breakdown_list<-c(bd31cmembership_mid_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd31cmembership_break_unique <-levels(as.factor(bd31cmembership[,2]))[-which(levels(as.factor(bd31cmembership[,2]))%in%unique(as.vector(replacement)))]
###add mid unique name to 0 

if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd31cmembership[which(bd31cmembership[,1]==levels(as.factor(bd31cmembership[,1]))[-1][which(colSums(overlap)==0)]),4:6]<-0
	}
}


TOM1_NEW_membership_para1<-bd31cmembership[,2]
TOM1_NEW_membership_para4<-bd31cmembership[,3]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31cmembership<-cbind(bd31cmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31cmembership<-cbind(bd31cmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31cmembership<-cbind(bd31cmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd31cmembership[,2])
bd31cmembership_bigger_breakdown_list<-bd31cmembership_bigger_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd31cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd31cmembership[which(bd31cmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd31cmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}


###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd31cmembership_bigger_breakdown_list<-c(bd31cmembership_bigger_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd31cmembership_break_unique <-c(bd31cmembership_break_unique ,levels(as.factor(bd31cmembership[,2]))[-which(levels(as.factor(bd31cmembership[,2]))%in%unique(as.vector(replacement)))])
###add bigger unique name to all at the end and bd31cmembership_bigger_unique
if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd31cmembership[which(bd31cmembership[,1]==levels(as.factor(bd31cmembership[,1]))[-1][which(colSums(overlap)==0)]),7:9]<-0
	}
}


##############################################################
bd31dmembership<-cbind(bd31dmembership0,bd31dmembershiploose,bd31dmembershipstrict)

TOM1_NEW_membership_para1<-bd31dmembership[,2]
TOM1_NEW_membership_para4<-bd31dmembership[,1]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31dmembership<-cbind(bd31dmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31dmembership<-cbind(bd31dmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31dmembership<-cbind(bd31dmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd31dmembership[,2])
bd31dmembership_mid_breakdown_list<-bd31dmembership_mid_unique <-bd31dmembership_break_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd31dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd31dmembership[which(bd31dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd31dmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}
max_number
###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd31dmembership_mid_breakdown_list<-c(bd31dmembership_mid_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd31dmembership_break_unique <-levels(as.factor(bd31dmembership[,2]))[-which(levels(as.factor(bd31dmembership[,2]))%in%unique(as.vector(replacement)))]
###add mid unique name to 0 

if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd31dmembership[which(bd31dmembership[,1]==levels(as.factor(bd31dmembership[,1]))[-1][which(colSums(overlap)==0)]),4:6]<-0
	}
}


TOM1_NEW_membership_para1<-bd31dmembership[,2]
TOM1_NEW_membership_para4<-bd31dmembership[,3]
moduleLabels11 = levels(as.factor(TOM1_NEW_membership_para1))[-1]
moduleLabels22 = levels(as.factor(TOM1_NEW_membership_para4))[-1]

# Numbers of modules
n1Mods = length(moduleLabels11)
n2Mods = length(moduleLabels22)
# Initialize tables of p-values and of the corresponding counts
pTbl = matrix(0, nrow = n1Mods, ncol = n2Mods);
# Execute all pairwaise comparisons
for (fmod in 1: n1Mods)
  for (cmod in 1: n2Mods)
  {
		genesets =list(a=which(TOM1_NEW_membership_para1 ==  moduleLabels11[fmod]),b=which(TOM1_NEW_membership_para4 ==  moduleLabels22[cmod]))
pTbl[fmod, cmod]<-MSET(genesets, n=total, lower.tail=FALSE)$p.value	
}
overlap = matrix(0, nrow = n1Mods, ncol = n2Mods);
overlap[which(pTbl <pcutoff)]<-1
###add names to 4,5,6
replacement<-matrix(0,dim(overlap)[2],4)
for (i in 1:dim(overlap)[2]){	
	if(length(which(overlap[,i]==1))>0){
	replacement[i,1:length(which(overlap[,i]==1))]<-which(overlap[,i]==1)
}
}
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
			replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31dmembership<-cbind(bd31dmembership, TOM1_NEW_membership_para4_rename)
replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>1){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[2]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31dmembership<-cbind(bd31dmembership, TOM1_NEW_membership_para4_rename)

replacement4<-rep(0,dim(replacement)[1])
for (i in 1:dim(replacement)[1]){
	if (length(which(unique(replacement[i,])>0))>2){
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[3]			
	}else{
		replacement4[i]<-unique(replacement[i,])[which(unique(replacement[i,])>0)]	[1]
		}
}
TOM1_NEW_membership_para4_rename<-
plyr::mapvalues(TOM1_NEW_membership_para4, from = c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1])), to = c(replacement4))
bd31dmembership<-cbind(bd31dmembership, TOM1_NEW_membership_para4_rename)

max_number<-max(bd31dmembership[,2])
bd31dmembership_bigger_breakdown_list<-bd31dmembership_bigger_unique <-NULL

###change duplicated name in mid to different 
if(length(which(rowSums(overlap)>1))>0){
for (i in 1:length(which(rowSums(overlap)>1))){
	for(j in 1:length(which(overlap[which(rowSums(overlap)>1)[i],]>0))){		
		print(max_number +1)
		print(c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]])
		print(which(bd31dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]))
bd31dmembership[which(bd31dmembership[,1]==c(as.numeric(levels(as.factor(TOM1_NEW_membership_para4))[-1]))[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j]]),4:6]<-max_number +1
print(levels(as.factor(bd31dmembership[,4])))
	replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],which(replacement[which(overlap[which(rowSums(overlap)>1)[i],]>0)[j],]==0)[1]]<-max_number +1	
		max_number= max_number+1
				
	}
}	
}


###break_breakdownlist
for (i in 1:dim(replacement)[1]){
if (length(unique(replacement[i,which(replacement[i,]>0)]))>1){
bd31dmembership_bigger_breakdown_list<-c(bd31dmembership_bigger_breakdown_list,toString(unique(replacement[i,which(replacement[i,]>0)])))
}
}
###breakunique
bd31dmembership_break_unique <-c(bd31dmembership_break_unique ,levels(as.factor(bd31dmembership[,2]))[-which(levels(as.factor(bd31dmembership[,2]))%in%unique(as.vector(replacement)))])
###add bigger unique name to all at the end and bd31dmembership_bigger_unique
if(length(which(colSums(overlap)==0))>0){
	for (i in 1:length(which(colSums(overlap)==0))){
	bd31dmembership[which(bd31dmembership[,1]==levels(as.factor(bd31dmembership[,1]))[-1][which(colSums(overlap)==0)]),7:9]<-0
	}
}


##generate reference network bigmatrix
##check this
list1_label=c(2,4,7)
list2_label=c(2,4,7)
list3_label=c(2,4,7)
list4_label=c(2,4,7)

list_i=0
list_all<-matrix(0,81,4)
for ( i1 in list1_label){
	for ( j1 in list2_label){
		for ( k1 in list3_label){
			for ( l1 in list4_label){
list_i<-list_i+1
list_all[list_i,]<-c(i1,j1,k1,l1)
}
}
}
}
print(list_all)
res<-function(input) {
#input<-list_all[1,]
bd21cmembership0<-c(bd21cmembership[,c(2,4:9)])
bd31cmembership0<-c(bd31cmembership[,c(2,4:9)])
bd21dmembership0<-c(bd21dmembership[,c(2,4:9)])
bd31dmembership0<-c(bd31dmembership[,c(2,4:9)])
threshold4<-threshold3<-threshold2<-threshold1<-0.05/length(levels(as.factor(bd21cmembership0)))/length(levels(as.factor(bd21dmembership0)))/length(levels(as.factor(bd31cmembership0)))/length(levels(as.factor(bd31dmembership0)))
a4=a3_1=a3_2=a3_3=a3_4=a2_1=a2_2=a2_3=a2_4=a2_5=a2_6=1
total=16207
list4<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])*length(levels(as.factor(bd21dmembership0))[-1])*length(levels(as.factor(bd31cmembership0))[-1])*length(levels(as.factor(bd31dmembership0))[-1])))
list3_1<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])*length(levels(as.factor(bd21dmembership0))[-1])*length(levels(as.factor(bd31cmembership0))[-1])))
list3_2<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])*length(levels(as.factor(bd21dmembership0))[-1])*length(levels(as.factor(bd31dmembership0))[-1])))
list3_3<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])*length(levels(as.factor(bd31cmembership0))[-1])*length(levels(as.factor(bd31dmembership0))[-1])))
list3_4<-matrix(0,1,(length(levels(as.factor(bd21dmembership0))[-1])*length(levels(as.factor(bd31cmembership0))[-1])*length(levels(as.factor(bd31dmembership0))[-1])))
list2_1<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])*length(levels(as.factor(bd21dmembership0))[-1])))
list2_2<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])*length(levels(as.factor(bd31cmembership0))[-1])))
list2_3<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])*length(levels(as.factor(bd31dmembership0))[-1])))
list2_4<-matrix(0,1,(length(levels(as.factor(bd21dmembership0))[-1])*length(levels(as.factor(bd31cmembership0))[-1])))
list2_5<-matrix(0,1,(length(levels(as.factor(bd21dmembership0))[-1])*length(levels(as.factor(bd31dmembership0))[-1])))
list2_6<-matrix(0,1,(length(levels(as.factor(bd31cmembership0))[-1])*length(levels(as.factor(bd31dmembership0))[-1])))
list1_1<-matrix(0,1,(length(levels(as.factor(bd21cmembership0))[-1])))
list1_2<-matrix(0,1,(length(levels(as.factor(bd21dmembership0))[-1])))
list1_3<-matrix(0,1,(length(levels(as.factor(bd31cmembership0))[-1])))
list1_4<-matrix(0,1,(length(levels(as.factor(bd31dmembership0))[-1])))
colnames(list4)<-rep(0,dim(list4)[2])
colnames(list3_1)<-rep(0,dim(list3_1)[2])
colnames(list3_2)<-rep(0,dim(list3_2)[2])
colnames(list3_3)<-rep(0,dim(list3_3)[2])
colnames(list3_4)<-rep(0,dim(list3_4)[2])
colnames(list2_1)<-rep(0,dim(list2_1)[2])
colnames(list2_2)<-rep(0,dim(list2_2)[2])
colnames(list2_3)<-rep(0,dim(list2_3)[2])
colnames(list2_4)<-rep(0,dim(list2_4)[2])
colnames(list2_5)<-rep(0,dim(list2_5)[2])
colnames(list2_6)<-rep(0,dim(list2_6)[2])
for (i in levels(as.factor(bd21cmembership0))[-1]){
	for (j in levels(as.factor(bd21dmembership0))[-1]){
		for (l in levels(as.factor(bd31cmembership0))[-1]){
			for (k in levels(as.factor(bd31dmembership0))[-1]){
				#i=j=l=k=1
if (paste(i,j,l,k)%in%colnames(list4)=="FALSE"){
		colnames(list4)[a4]<-paste(i,j,l,k)
		a4<-a4+1
		}
if (paste(i,j,l,"-")%in%colnames(list3_1)=="FALSE"){
		colnames(list3_1)[a3_1]<-paste(i,j,l,"-")
		a3_1<-a3_1+1
		}
if (paste(i,j,"-",k)%in%colnames(list3_2)=="FALSE"){
		colnames(list3_2)[a3_2]<-paste(i,j,"-",k)
		a3_2<-a3_2+1
		}
if (paste(i,"-",l,k)%in%colnames(list3_3)=="FALSE"){
		colnames(list3_3)[a3_3]<-paste(i,"-",l,k)
		a3_3<-a3_3+1
		}
if (paste("-",j,l,k)%in%colnames(list3_4)=="FALSE"){
		colnames(list3_4)[a3_4]<-paste("-",j,l,k)
		a3_4<-a3_4+1
		}
if (paste(i,j,"-","-")%in%colnames(list2_1)=="FALSE"){
		colnames(list2_1)[a2_1]<-paste(i,j,"-","-")
		a2_1 <-a2_1 +1
		}
if (paste(i,"-",l,"-")%in%colnames(list2_2)=="FALSE"){
		colnames(list2_2)[a2_2]<-paste(i,"-",l,"-")
		a2_2 <-a2_2 +1
		}
if (paste(i,"-","-",k)%in%colnames(list2_3)=="FALSE"){
		colnames(list2_3)[a2_3]<-paste(i,"-","-",k)
		a2_3 <-a2_3 +1
		}		
if (paste("-",j,l,"-")%in%colnames(list2_4)=="FALSE"){
		colnames(list2_4)[a2_4]<-paste("-",j,l,"-")
		a2_4 <-a2_4 +1
		}
if (paste("-",j,"-",k)%in%colnames(list2_5)=="FALSE"){
		colnames(list2_5)[a2_5]<-paste("-",j,"-",k)
		a2_5 <-a2_5 +1
		}	
if (paste("-","-",l,k)%in%colnames(list2_6)=="FALSE"){
		colnames(list2_6)[a2_6]<-paste("-","-",l,k)
		a2_6 <-a2_6 +1
		}	
}
}
}	
}


colnames(list1_1)<-paste(levels(as.factor(bd21cmembership0))[-1],"-","-","-")
colnames(list1_2)<-paste("-",levels(as.factor(bd21dmembership0))[-1],"-","-")
colnames(list1_3)<-paste("-","-",levels(as.factor(bd31cmembership0))[-1],"-")
colnames(list1_4)<-paste("-","-","-",levels(as.factor(bd31dmembership0))[-1])
list4_number<-list4
list1_1_number<-list1_1
list1_2_number<-list1_2
list1_3_number<-list1_3
list1_4_number<-list1_4
list2_1_number<-list2_1
list2_2_number<-list2_2
list2_3_number<-list2_3
list2_4_number<-list2_4
list2_5_number<-list2_5
list2_6_number<-list2_6
list3_1_number<-list3_1
list3_2_number<-list3_2
list3_3_number<-list3_3
list3_4_number<-list3_4

bd21cmembership0<-bd21cmembership[,input[1]]
bd21dmembership0<-bd21dmembership[,input[2]]
bd31cmembership0<-bd31cmembership[,input[3]]
bd31dmembership0<-bd31dmembership[,input[4]]
ain81<-1
for (i in levels(as.factor(bd21cmembership0))[-1]){
	for (j in levels(as.factor(bd21dmembership0))[-1]){
	for (l in levels(as.factor(bd31cmembership0))[-1]){
			for (k in levels(as.factor(bd31dmembership0))[-1]){
#i=j=l=7
fit<-fit1<-fit2<-fit3<-fit4<-fit1_1<-fit1_2<-fit1_3<-fit1_4<-fit1_5<-fit1_6<-NULL
iextend=i
jextend=j 	
lextend=l
kextend=k	
genesets<-list(a=colnames(matrixbd)[which(bd21cmembership0 ==i)],b=colnames(matrixbd)[which(bd21dmembership0 ==j)],c=colnames(matrixbd)[which(bd31cmembership0 ==l)],d=colnames(matrixbd)[which(bd31dmembership0 ==k)])
fit=MSET(genesets, n=total, lower.tail=FALSE)
if ((fit$p.value < threshold4)){
 extendlist	<-apply(as.matrix(expand.grid(as.character(iextend),as.character(jextend),as.character(lextend),as.character(kextend))), 1, paste, collapse=" ")
list4[ain81,which(colnames(list4)%in% extendlist)]<-ifelse(fit$p.value==0,10,fit$p.value)
  	list4_number[ain81,which(colnames(list4)%in% extendlist)]<-length(fit$intersects)
}
genesets1 <-list(a=colnames(matrixbd)[which(bd21cmembership0 ==i)],b=colnames(matrixbd)[which(bd21dmembership0 ==j)],c=colnames(matrixbd)[which(bd31cmembership0 ==l)])
fit1=MSET(genesets1, n=total, lower.tail=FALSE)
genesets2<-list(a=colnames(matrixbd)[which(bd21cmembership0 ==i)],b=colnames(matrixbd)[which(bd21dmembership0 ==j)],c=colnames(matrixbd)[which(bd31dmembership0 ==k)])
fit2=MSET(genesets2, n=total, lower.tail=FALSE)
genesets3<-list(a=colnames(matrixbd)[which(bd21cmembership0 ==i)],b=colnames(matrixbd)[which(bd31cmembership0 ==l)],c=colnames(matrixbd)[which(bd31dmembership0 ==k)])
fit3=MSET(genesets3, n=total, lower.tail=FALSE)
genesets4<-list(a=colnames(matrixbd)[which(bd21dmembership0 ==j)],b=colnames(matrixbd)[which(bd31cmembership0 ==l)],c=colnames(matrixbd)[which(bd31dmembership0 ==k)])
fit4=MSET(genesets4, n=total, lower.tail=FALSE)
if (((fit1$p.value < threshold4))){
 extendlist	<-apply(as.matrix(expand.grid(as.character(iextend),as.character(jextend),as.character(lextend),"-")), 1, paste, collapse=" ")
  	list3_1[ain81,which(colnames(list3_1)%in% extendlist)]<-ifelse(fit1$p.value==0,10,fit1$p.value)
  	 	list3_1_number[ain81,which(colnames(list3_1)%in% extendlist)]<-length(fit1$intersects)
}
if(((fit2$p.value < threshold4))){
 	 extendlist	<-apply(as.matrix(expand.grid(as.character(iextend),as.character(jextend),"-",as.character(kextend))), 1, paste, collapse=" ")	
 	list3_2[ain81,which(colnames(list3_2)%in% extendlist)]<-ifelse(fit2$p.value==0,10,fit2$p.value)
 	  	 	list3_2_number[ain81,which(colnames(list3_2)%in% extendlist)]<-length(fit2$intersects)
}
if (((fit3$p.value < threshold4))){
	 	 extendlist	<-apply(as.matrix(expand.grid(as.character(iextend),"-",as.character(lextend),as.character(kextend))), 1, paste, collapse=" ")
 	list3_3[ain81,which(colnames(list3_3)%in% extendlist)]<-ifelse(fit3$p.value==0,10,fit3$p.value)
 	list3_3_number[ain81,which(colnames(list3_3)%in% extendlist)]<-length(fit3$intersects)
}
if (((fit4$p.value < threshold4))){
	 	 extendlist	<-apply(as.matrix(expand.grid("-",as.character(jextend),as.character(lextend),as.character(kextend))), 1, paste, collapse=" ")
 	list3_4[ain81,which(colnames(list3_4)%in% extendlist)]<-ifelse(fit4$p.value==0,10,fit4$p.value)
 	list3_4_number[ain81,which(colnames(list3_4)%in% extendlist)]<-length(fit4$intersects)
}
genesets1_1<-list(a=colnames(matrixbd)[which(bd21cmembership0 ==i)],b=colnames(matrixbd)[which(bd21dmembership0 ==j)])
fit1_1=MSET(genesets1_1, n=total, lower.tail=FALSE)
genesets1_2<-list(a=colnames(matrixbd)[which(bd21cmembership0 ==i)],b=colnames(matrixbd)[which(bd31cmembership0 ==l)])
fit1_2=MSET(genesets1_2, n=total, lower.tail=FALSE)
genesets1_3<-list(a=colnames(matrixbd)[which(bd21cmembership0 ==i)],b=colnames(matrixbd)[which(bd31dmembership0 ==k)])
fit1_3=MSET(genesets1_3, n=total, lower.tail=FALSE)
genesets1_4<-list(a=colnames(matrixbd)[which(bd21dmembership0 ==j)],b=colnames(matrixbd)[which(bd31cmembership0 ==l)])
fit1_4=MSET(genesets1_4, n=total, lower.tail=FALSE)
genesets1_5<-list(a=colnames(matrixbd)[which(bd21dmembership0 ==j)],b=colnames(matrixbd)[which(bd31dmembership0 ==k)])
fit1_5=MSET(genesets1_5, n=total, lower.tail=FALSE)
genesets1_6<-list(a=colnames(matrixbd)[which(bd31cmembership0 ==l)],b=colnames(matrixbd)[which(bd31dmembership0 ==k)])
fit1_6=MSET(genesets1_6, n=total, lower.tail=FALSE)
if (fit1_1 $p.value<threshold2){
	 	 extendlist	<-apply(as.matrix(expand.grid(as.character(iextend),as.character(jextend),"-","-")), 1, paste, collapse=" ")
 	list2_1[ain81,which(colnames(list2_1)%in% extendlist)]<-ifelse(fit1_1$p.value==0,10,fit1_1$p.value)
 	list2_1_number[ain81,which(colnames(list2_1)%in% extendlist)]<-length(fit1_1$intersects)
}
if (fit1_2 $p.value<threshold2){
	 	 extendlist	<-apply(as.matrix(expand.grid(as.character(iextend),"-",as.character(lextend),"-")), 1, paste, collapse=" ")
 	list2_2[ain81,which(colnames(list2_2)%in% extendlist)]<-ifelse(fit1_2$p.value==0,10,fit1_2$p.value)
 	 	list2_2_number[ain81,which(colnames(list2_2)%in% extendlist)]<-length(fit1_2$intersects) 	
 	}
if (fit1_3 $p.value<threshold2){
	 	 extendlist	<-apply(as.matrix(expand.grid(as.character(iextend),"-","-",as.character(kextend))), 1, paste, collapse=" ")
 	list2_3[ain81,which(colnames(list2_3)%in% extendlist)]<-ifelse(fit1_3$p.value==0,10,fit1_3$p.value)
 	 	list2_3_number[ain81,which(colnames(list2_3)%in% extendlist)]<-length(fit1_3$intersects)
}
if (fit1_4 $p.value<threshold2){
	 	 extendlist	<-apply(as.matrix(expand.grid("-",as.character(jextend),as.character(lextend),"-")), 1, paste, collapse=" ")
 	list2_4[ain81,which(colnames(list2_4)%in% extendlist)]<-ifelse(fit1_4$p.value==0,10,fit1_4$p.value)
 	 	list2_4_number[ain81,which(colnames(list2_4)%in% extendlist)]<-length(fit1_4$intersects)
 	
}
if (fit1_5 $p.value<threshold2){
	 	 extendlist	<-apply(as.matrix(expand.grid("-",as.character(jextend),"-",as.character(kextend))), 1, paste, collapse=" ")
 	list2_5[ain81,which(colnames(list2_5)%in% extendlist)]<-ifelse(fit1_5$p.value==0,10,fit1_5$p.value)
 	 	list2_5_number[ain81,which(colnames(list2_5)%in% extendlist)]<-length(fit1_5$intersects)
 		
}
if (fit1_6 $p.value<threshold2){
	 	 extendlist	<-apply(as.matrix(expand.grid("-","-",as.character(lextend),as.character(kextend))), 1, paste, collapse=" ")
 	list2_6[ain81,which(colnames(list2_6)%in% extendlist)]<-ifelse(fit1_6$p.value==0,10,fit1_6$p.value)
 	 	list2_6_number[ain81,which(colnames(list2_6)%in% extendlist)]<-length(fit1_6$intersects)
}

}
}	
}
}	
output<-rbind(c(list1_1[1,], list1_2[1,], list1_3[1,], list1_4[1,], list2_1[1,], list2_2[1,], list2_3[1,], list2_4[1,], list2_5[1,], list2_6[1,], list3_1[1,], list3_2[1,], list3_3[1,], list3_4[1,], list4[1,]),c(list1_1_number[1,], list1_2_number[1,], list1_3_number[1,], list1_4_number[1,], list2_1_number[1,], list2_2_number[1,], list2_3_number[1,], list2_4_number[1,], list2_5_number[1,], list2_6_number[1,], list3_1_number[1,], list3_2_number[1,], list3_3_number[1,], list3_4_number[1,], list4_number[1,]))
#output<-c(test1,test2,test3,test4, test5, test6, test7, test8,levels(as.factor(bd21cmembership0))[-1],levels(as.factor(bd21dmembership0))[-1],levels(as.factor(bd31cmembership0))[-1],levels(as.factor(bd31dmembership0))[-1])
return(output)
}

args = commandArgs(trailingOnly=TRUE)
print(args)
print(args[1])
it<-as.numeric(args[1])
finalmatirx<-res(list_all[it,])
finalmatirx_1 <-finalmatirx[,which(colSums(finalmatirx)>0)]
head(finalmatirx_1)

write.csv(finalmatirx_1[1,], file=paste("/home/gridsan/jyun/network/stable_network_finalmatirx_12_15",args[1],".csv",sep=""))
write.csv(finalmatirx_1[2,], file=paste("/home/gridsan/jyun/network/stable_network_finalmatirx_number_12_15",args[1],".csv",sep=""))

