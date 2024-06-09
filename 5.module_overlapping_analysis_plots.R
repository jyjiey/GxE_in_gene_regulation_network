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
overlap[which(pTbl <10e-13)]<-1
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
overlap[which(pTbl <10e-13)]<-1
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
overlap[which(pTbl <10e-13)]<-1
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
overlap[which(pTbl <10e-13)]<-1
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
overlap[which(pTbl <10e-13)]<-1
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
overlap[which(pTbl <10e-13)]<-1
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
overlap[which(pTbl <10e-13)]<-1
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
overlap[which(pTbl <10e-13)]<-1
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
 #low cutoff=35, high cutoff= 80 low p=0.00000001, high #p=0.00000000000000001

#p=0.000000000001
p=p3= p2=1e-12

#install.packages("qlcMatrix",repos="https://cran.r-project.org")

library(qlcMatrix)
##why different? 
for (k in 1:81){
test2<-read.csv(paste("/home/gridsan/jyun/network/stable_network_finalmatirx_12_15",k,".csv",sep=""), header=T,row.names=1)
test1<-read.csv(paste("/home/gridsan/jyun/network/stable_network_finalmatirx_number_12_15",k,".csv",sep=""), header=T,row.names=1)

#test2<-read.csv(paste("/home/gridsan/jyun/network/stable_network_finalmatirx_10_15",k,".csv",sep=""), header=T,row.names=1)
#test1<-read.csv(paste("/home/gridsan/jyun/network/stable_network_finalmatirx_number_10_15",k,".csv",sep=""), header=T,row.names=1)
singleton<-matrix(1,2,length(c(colnames(list1_1),colnames(list1_2),colnames(list1_3),colnames(list1_4))))
tem<-cbind(test1,test2)
for (i in 1:dim(cbind(test1,test2))[1]){
	if (length(which(strsplit (rownames(cbind(test1,test2)[i,]), split = " ")[[1]]=="-"))==1){
		tem[i,2]= ifelse((tem[i,2]<p3)|(tem[i,2]==10),tem[i,2],1)		
	}	
}
test2 <-tem[,2]
for (i in 1:dim(cbind(test1,test2))[1]){
	if (length(which(strsplit (rownames(cbind(test1,test2)[i,]), split = " ")[[1]]=="-"))==0){
		tem[i,2]= ifelse((tem[i,2]<p)|(tem[i,2]==10),tem[i,2],1)		
	}	
}
test2 <-tem[,2]
for (i in 1:dim(cbind(test1,test2))[1]){
	if (length(which(strsplit (rownames(cbind(test1,test2)[i,]), split = " ")[[1]]=="-"))==2){
		tem[i,2]= ifelse((tem[i,2]<p2)|(tem[i,2]==10),tem[i,2],1)		
	}	
}
test2 <-tem[,2]

colnames(singleton)<-c(colnames(list1_1),colnames(list1_2),colnames(list1_3),colnames(list1_4))
finalmatirx_1<-t(cbind(test1,test2)[rev(which((((test2<p2)|(test2>1))&(test1>1))|test2==10)),])
finalmatirx_1<-cbind(finalmatirx_1, singleton)


dim(finalmatirx_1)

bd21cmembership_all_breakdown_list <-c(bd21cmembership_mid_breakdown_list, bd21cmembership_bigger_breakdown_list)
duplicatedlist<-c(as.numeric(strsplit(toString(bd21cmembership_mid_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd21cmembership_mid_breakdown_list), split = ",")[[1]]))],as.numeric(strsplit(toString(bd21cmembership_bigger_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd21cmembership_bigger_breakdown_list), split = ",")[[1]]))])
finalmatirx_1_fix<-finalmatirx_1

if(length(bd21cmembership_mid_breakdown_list)>0){
for (j in 1:length(bd21cmembership_mid_breakdown_list)){
	list<-as.numeric(strsplit(bd21cmembership_mid_breakdown_list[j], split = ",")[[1]])
	
	for (ii in 1:dim(finalmatirx_1_fix)[2]){
		if(length(which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii]))>0){
		i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	if (strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1]%in%list){


		newnames<-paste(list ,strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2], strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4])
		
	add<-as.vector(qlcMatrix::rowMax(finalmatirx_1 [,which(colnames(finalmatirx_1)%in%newnames)],ignore.zero =FALSE))
finalmatirx_1 <-cbind(add,finalmatirx_1)


	if(length(which(list %in% duplicatedlist))>0){

		newnames<-newnames[-which(newnames%in%paste(list[which(list%in% duplicatedlist)],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4]))]
		duplicatedlist <-duplicatedlist[-which(duplicatedlist%in%list)]
	}
	i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	colnames(finalmatirx_1)[1]<-paste(list[which(list %in%bd21cmembership[,4])] ,strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2], strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4])
if(length(which(colnames(finalmatirx_1)%in%newnames))>1){	
finalmatirx_1 <- finalmatirx_1[,-which(colnames(finalmatirx_1)%in%newnames)[-1]]

	}
}
}
}
}
}
dim(finalmatirx_1)

dim(finalmatirx_1)

bd21dmembership_all_breakdown_list <-c(bd21dmembership_mid_breakdown_list, bd21dmembership_bigger_breakdown_list)
duplicatedlist<-c(as.numeric(strsplit(toString(bd21dmembership_mid_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd21dmembership_mid_breakdown_list), split = ",")[[1]]))],as.numeric(strsplit(toString(bd21dmembership_bigger_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd21dmembership_bigger_breakdown_list), split = ",")[[1]]))])
finalmatirx_1_fix<-finalmatirx_1

if(length(bd21dmembership_mid_breakdown_list)>0){
for (j in 1:length(bd21dmembership_mid_breakdown_list)){
	list<-as.numeric(strsplit(bd21dmembership_mid_breakdown_list[j], split = ",")[[1]])
	
	for (ii in 1:dim(finalmatirx_1_fix)[2]){
		if(length(which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii]))>0){
		i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	if (strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2]%in%list){


		newnames<-paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1],list , strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4])
		
	add<-as.vector(qlcMatrix::rowMax(finalmatirx_1 [,which(colnames(finalmatirx_1)%in%newnames)],ignore.zero =FALSE))
finalmatirx_1 <-cbind(add,finalmatirx_1)


	if(length(which(list %in% duplicatedlist))>0){

		newnames<-newnames[-which(newnames%in%paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1],list[which(list%in% duplicatedlist)],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4]))]
		duplicatedlist <-duplicatedlist[-which(duplicatedlist%in%list)]
	}
	i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	colnames(finalmatirx_1)[1]<-paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1],list[which(list %in%bd21dmembership[,4])] , strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4])
if(length(which(colnames(finalmatirx_1)%in%newnames))>1){	
finalmatirx_1 <- finalmatirx_1[,-which(colnames(finalmatirx_1)%in%newnames)[-1]]

	}
}
}
}
}
}
dim(finalmatirx_1)


dim(finalmatirx_1)

bd31cmembership_all_breakdown_list <-c(bd31cmembership_mid_breakdown_list, bd31cmembership_bigger_breakdown_list)
duplicatedlist<-c(as.numeric(strsplit(toString(bd31cmembership_mid_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd31cmembership_mid_breakdown_list), split = ",")[[1]]))],as.numeric(strsplit(toString(bd31cmembership_bigger_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd31cmembership_bigger_breakdown_list), split = ",")[[1]]))])
finalmatirx_1_fix<-finalmatirx_1

if(length(bd31cmembership_mid_breakdown_list)>0){
for (j in 1:length(bd31cmembership_mid_breakdown_list)){
	list<-as.numeric(strsplit(bd31cmembership_mid_breakdown_list[j], split = ",")[[1]])
	
	for (ii in 1:dim(finalmatirx_1_fix)[2]){
		if(length(which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii]))>0){
		i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	if (strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3]%in%list){


		newnames<-paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1], strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2],list ,strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4])
		
	add<-as.vector(qlcMatrix::rowMax(finalmatirx_1 [,which(colnames(finalmatirx_1)%in%newnames)],ignore.zero =FALSE))
finalmatirx_1 <-cbind(add,finalmatirx_1)


	if(length(which(list %in% duplicatedlist))>0){

		newnames<-newnames[-which(newnames%in%paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2],list[which(list%in% duplicatedlist)],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4]))]
		duplicatedlist <-duplicatedlist[-which(duplicatedlist%in%list)]
	}
	i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	colnames(finalmatirx_1)[1]<-paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1], strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2],list[which(list %in%bd31cmembership[,4])] ,strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4])
if(length(which(colnames(finalmatirx_1)%in%newnames))>1){	
finalmatirx_1 <- finalmatirx_1[,-which(colnames(finalmatirx_1)%in%newnames)[-1]]

	}
}
}
}
}
}
dim(finalmatirx_1)
dim(finalmatirx_1)

bd31dmembership_all_breakdown_list <-c(bd31dmembership_mid_breakdown_list, bd31dmembership_bigger_breakdown_list)
duplicatedlist<-c(as.numeric(strsplit(toString(bd31dmembership_mid_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd31dmembership_mid_breakdown_list), split = ",")[[1]]))],as.numeric(strsplit(toString(bd31dmembership_bigger_breakdown_list), split = ",")[[1]])[duplicated(as.numeric(strsplit(toString(bd31dmembership_bigger_breakdown_list), split = ",")[[1]]))])
finalmatirx_1_fix<-finalmatirx_1

if(length(bd31dmembership_mid_breakdown_list)>0){
for (j in 1:length(bd31dmembership_mid_breakdown_list)){
	list<-as.numeric(strsplit(bd31dmembership_mid_breakdown_list[j], split = ",")[[1]])
	
	for (ii in 1:dim(finalmatirx_1_fix)[2]){
		if(length(which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii]))>0){
		i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	if (strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][4]%in%list){


		newnames<-paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1], strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],list )
		
	add<-as.vector(qlcMatrix::rowMax(finalmatirx_1 [,which(colnames(finalmatirx_1)%in%newnames)],ignore.zero =FALSE))
finalmatirx_1 <-cbind(add,finalmatirx_1)


	if(length(which(list %in% duplicatedlist))>0){

		newnames<-newnames[-which(newnames%in%paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],list[which(list%in% duplicatedlist)]))]
		duplicatedlist <-duplicatedlist[-which(duplicatedlist%in%list)]
	}
	i<-which(colnames(finalmatirx_1)==colnames(finalmatirx_1_fix)[ii])
	colnames(finalmatirx_1)[1]<-paste(strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][1], strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit (colnames(finalmatirx_1)[i], split = " ")[[1]][3],list[which(list %in%bd31dmembership[,4])] )
	print(list)
	print(list[which(list %in%bd31dmembership[,4])] )
	print(colnames(finalmatirx_1)[1])
if(length(which(colnames(finalmatirx_1)%in%newnames))>1){	
finalmatirx_1 <- finalmatirx_1[,-which(colnames(finalmatirx_1)%in%newnames)[-1]]

	}
}
}
}
}
}
dim(finalmatirx_1)

finalmatirx_1_tem<-finalmatirx_1
for (j in 1:dim(finalmatirx_1_tem)[2]){
	print(i)
	if(colnames(finalmatirx_1_tem)[j]%in%colnames(finalmatirx_1)){
		i=which(colnames(finalmatirx_1)==colnames(finalmatirx_1_tem)[j])

names<-c(which(colnames(finalmatirx_1)%in%c(paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]),paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],"-"),paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],"-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]),paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],"-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]),paste("-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]),paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],"-","-"),paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],"-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],"-"),paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],"-","-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]),paste("-", strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],"-"),paste("-", strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],"-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]),paste("-","-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]),paste(strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][1],"-","-","-"),paste("-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][2],"-","-"),paste("-","-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][3],"-"),paste("-","-","-",strsplit(colnames(finalmatirx_1)[i], split = " ")[[1]][4]))))	
#print(colnames(finalmatirx_1)[names])

if(length(names)>1){
	print("--")
	print(names)
	print(i)
	newnames<-names[-which(names==i)]
	#print(dim(finalmatirx_1))
		#print(colnames(finalmatirx_1)[names])
finalmatirx_1<-finalmatirx_1[,-newnames]
		#print(dim(finalmatirx_1))

}


}
}
#print(t(finalmatirx_1))
write.csv(finalmatirx_1[1,], file=paste("/home/gridsan/jyun/network/stable_network_finalmatirx_12_15real",k,".csv",sep=""))
write.csv(finalmatirx_1[2,], file=paste("/home/gridsan/jyun/network/stable_network_finalmatirx_number_12_15real",k,".csv",sep=""))
}

###together
finalmatirx_1<-cbind(list1_1, list1_2, list1_3, list1_4, list2_1, list2_2, list2_3, list2_4, list2_5, list2_6, list3_1, list3_2, list3_3, list3_4, list4)
finalmatirx_1_number<-finalmatirx_1
setwd("/home/gridsan/jyun/network/")
file.names <- dir("/home/gridsan/jyun/network/",pattern ="stable_network_finalmatirx_number_12_15real")
file.names2 <- dir("/home/gridsan/jyun/network/",pattern ="stable_network_finalmatirx_12_15real")

name<-matrix(0,1,81)
finalmatirx_1_tem<-finalmatirx_1_number_tem<-NULL
for (i in 1:length(file.names)) {	
matrix <- read.table(file.names[i],header= TRUE,sep=",")
matrix_number <- read.table(file.names2[i],header= TRUE,sep=",")
finalmatirx_1<-rep(0,length(finalmatirx_1))
finalmatirx_1[match(matrix$X,colnames(cbind(list1_1, list1_2, list1_3, list1_4, list2_1, list2_2, list2_3, list2_4, list2_5, list2_6, list3_1, list3_2, list3_3, list3_4, list4)))]<-matrix$x
finalmatirx_1_tem<-cbind(finalmatirx_1_tem, finalmatirx_1)
finalmatirx_1_number<-rep(0, length(finalmatirx_1_number))
finalmatirx_1_number[match(matrix_number $X,colnames(cbind(list1_1, list1_2, list1_3, list1_4, list2_1, list2_2, list2_3, list2_4, list2_5, list2_6, list3_1, list3_2, list3_3, list3_4, list4)))]<-matrix_number $x
finalmatirx_1_number_tem<-cbind(finalmatirx_1_number_tem, finalmatirx_1_number)
print(dim(matrix_number))
}
#finalmatirx_1[i,match(matrix$X,colnames(finalmatirx_1))]<-matrix$x
#finalmatirx_1_number[i,match(matrix_number $X,colnames(finalmatirx_1_number))]<-matrix_number $x

finalmatirx_1_tem_0<-t(finalmatirx_1_tem[which(rowSums(finalmatirx_1_tem)>0),])
finalmatirx_1_number_0<-t(finalmatirx_1_number_tem[which(rowSums(finalmatirx_1_number_tem)>0),])


colnames(finalmatirx_1_tem_0)<-colnames(finalmatirx_1_number_0)<-colnames(cbind(list1_1, list1_2, list1_3, list1_4, list2_1, list2_2, list2_3, list2_4, list2_5, list2_6, list3_1, list3_2, list3_3, list3_4, list4))[which(rowSums(finalmatirx_1_tem)>0)]
colSums(ifelse(finalmatirx_1_number_0>0,1,0))
##check <2 genes in 4 overlapping???
##category
#finalmatirx_1_number_keep1<-finalmatirx_1_number_0[ ,which(colSums(ifelse(finalmatirx_1_number_0>0,1,0))>cutoff)]
#colnames(finalmatirx_1_number_keep1)<-colnames(finalmatirx_1_tem_0)[which(colSums(ifelse(finalmatirx_1_number_0>0,1,0))>cutoff)]

###81
cutoff=80
finalmatirx_1_number_keep1<-finalmatirx_1_number_0[ ,which(colSums(ifelse(finalmatirx_1_number_0>0,1,0))>cutoff)]
colnames(finalmatirx_1_number_keep1)<-colnames(finalmatirx_1_tem_0)[which(colSums(ifelse(finalmatirx_1_number_0>0,1,0))>cutoff)]


finalmatirx_1_keep1<-finalmatirx_1_tem_0[ ,which(colSums(ifelse(finalmatirx_1_number_0>0,1,0))>cutoff)]
colnames(finalmatirx_1_keep1)<-colnames(finalmatirx_1_tem_0)[which(colSums(ifelse(finalmatirx_1_number_0>0,1,0))>cutoff)]


finalmatirx_1_number_keep<-finalmatirx_1_number_keep1
write.csv(finalmatirx_1_number_keep, file=paste("/home/gridsan/jyun/network/finalmatirx_1_number_keep_all.csv",sep=""))
write.csv(finalmatirx_1_keep1, file=paste("/home/gridsan/jyun/network/finalmatirx_1_pvalue_keep_all.csv",sep=""))
#addon<-rep(9,81)
#finalmatirx_1_number_keep<-cbind(finalmatirx_1_number_keep,rep(9,81))
#colnames(finalmatirx_1_number_keep)[length(colnames(finalmatirx_1_number_keep))]<-"17 16 15 19"


total_number<-c((length(levels(as.factor(bd21cmembership[,4]))[-1]))+(length(levels(as.factor(bd21dmembership[,4]))[-1]))+(length(levels(as.factor(bd31cmembership[,4]))[-1]))+(length(levels(as.factor(bd31dmembership[,4]))[-1])))+4
adjacency_prepare<-colnames(finalmatirx_1_number_keep)
	adjacency<-matrix(0, total_number, total_number)
	rownames(adjacency)<-colnames(adjacency)<-c(paste0("Bd21c-",c(levels(as.factor(bd21cmembership[,4]))[-1],"-")),paste0("Bd21d-",c(levels(as.factor(bd21dmembership[,4]))[-1],"-")),paste0("Bd31c-",c(levels(as.factor(bd31cmembership[,4]))[-1],"-")),paste0("Bd31d-",c(levels(as.factor(bd31dmembership[,4]))[-1],"-")))
	adjacency_prepare<-matrix(0,dim(finalmatirx_1_number_keep)[2],4)
	for (i in 1:dim(finalmatirx_1_number_keep)[2]){
	adjacency_prepare[i,]<-c(paste0("Bd21c-",strsplit(colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][1]),paste0("Bd21d-",strsplit(colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][2]),paste0("Bd31c-",strsplit(colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][3]),paste0("Bd31d-",strsplit(colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][4]))
	adjacency[which(colnames(adjacency)%in% adjacency_prepare[i,]),which(colnames(adjacency)%in% adjacency_prepare[i,])]<-1
}

adjacency<-adjacency[-which(colnames(adjacency)%in%c("Bd21c--","Bd21d--","Bd31c--","Bd31d--")),-which(colnames(adjacency)%in%c("Bd21c--","Bd21d--","Bd31c--","Bd31d--"))]
##network of strict overlapping?

######

 oldnames<-colnames(adjacency)
size1<-plyr::count(bd21cmembership[,4])[-1,2]
size2<-plyr::count(bd21dmembership[,4])[-1,2]
size3<-plyr::count(bd31cmembership[,4])[-1,2]
size4<-plyr::count(bd31dmembership[,4])[-1,2]
size<-c(size1,size2,size3,size4)
##addunique size

names(size)<-rownames(adjacency)

vertexcolor=c(rep("blue",(length(levels(as.factor(bd21cmembership[,4]))[-1]))),rep("red",(length(levels(as.factor(bd21dmembership[,4]))[-1]))),rep("purple",(length(levels(as.factor(bd31cmembership[,4]))[-1]))),rep("orange",(length(levels(as.factor(bd31dmembership[,4]))[-1]))))
names(vertexcolor)<-rownames(adjacency)


bd21cmembership_mid_breakdown_list_first <-bd21dmembership_mid_breakdown_list_first<-bd31cmembership_mid_breakdown_list_first<-bd31dmembership_mid_breakdown_list_first<-NULL

bd21cmembership_mid_breakdown_list_first <-as.numeric(strsplit(toString(bd21cmembership_mid_breakdown_list),split = ",")[[1]])[which(as.numeric(strsplit(toString(bd21cmembership_mid_breakdown_list),split = ",")[[1]])%in%bd21cmembership[,4])]
bd21dmembership_mid_breakdown_list_first <-as.numeric(strsplit(toString(bd21dmembership_mid_breakdown_list),split = ",")[[1]])[which(as.numeric(strsplit(toString(bd21dmembership_mid_breakdown_list),split = ",")[[1]])%in%bd21dmembership[,4])]

bd31cmembership_mid_breakdown_list_first <-as.numeric(strsplit(toString(bd31cmembership_mid_breakdown_list),split = ",")[[1]])[which(as.numeric(strsplit(toString(bd31cmembership_mid_breakdown_list),split = ",")[[1]])%in%bd31cmembership[,4])]

bd31dmembership_mid_breakdown_list_first <-as.numeric(strsplit(toString(bd31dmembership_mid_breakdown_list),split = ",")[[1]])[which(as.numeric(strsplit(toString(bd31dmembership_mid_breakdown_list),split = ",")[[1]])%in%bd31dmembership[,4])]


adjacency_prepare<-adjacency_prepare[-4,]
finalmatirx_1_number_keep_mark<-finalmatirx_1_number_keep
finalmatirx_1_number_keep_mark<-finalmatirx_1_number_keep_mark[,-4]
ID<-GROUP<-group_color <-group_color_FULL<-transparentness <-NULL
ID<-as.vector(t(adjacency_prepare))[-which(as.vector(t(adjacency_prepare))%in%c("Bd31c--","Bd31d--","Bd21c--", "Bd21d--"))]
for (i in 1:dim(adjacency_prepare)[1]){
GROUP<-c(GROUP,rep(paste0(c(4-length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd31c--","Bd31d--","Bd21c--", "Bd21d--"))])),"-",i),c(4-length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd31c--","Bd31d--","Bd21c--", "Bd21d--"))]))))
if (length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd31c--","Bd31d--","Bd21c--", "Bd21d--"))])==0){
group_color<-c(group_color,"#26C6DA")
}
if (length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd31c--","Bd31d--","Bd21c--", "Bd21d--"))])==1){
group_color<-c(group_color,"#EC407A")
}
if (length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd31c--","Bd31d--","Bd21c--", "Bd21d--"))])==2){
group_color<-c(group_color,"#4CAF50")
if (length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd21c--", "Bd21d--"))])==2){
group_color[i]<-"#FFFF00"
}
if (length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd31c--","Bd31d--"))])==2){
group_color[i]<-"#FFFF00"
}

}

if (length(adjacency_prepare[i,][which(adjacency_prepare[i,]%in%c("Bd31c--","Bd31d--","Bd21c--", "Bd21d--"))])==3){
group_color<-c(group_color,"#FFAB91")
}
transparentness<-colMax(finalmatirx_1_number_keep_mark[,i])
transparentness<-ifelse(transparentness>99,99, transparentness)
transparentness<-ifelse(transparentness<10,10, transparentness)
group_color_FULL<-c(group_color_FULL,paste0(group_color[i], transparentness))
}
members<-data_frame(id= ID,group = GROUP)
group_ids <- lapply(members %>% split(.$group), function(grp) { grp$id })
class(group_ids[1])

library(igraph)

Newbigmatrix1<-adjacency
Newbigmatrix1 <-adjacency[-which(rowSums(adjacency)==0),-which(rowSums(adjacency)==0)]
diag(Newbigmatrix1) <- 0
#
shape<-ifelse(gxe==1,"square","circle")
shape<-shape[-which(rowSums(adjacency)==0)]
STABLEgraph1 <-graph_from_adjacency_matrix(as.matrix(Newbigmatrix1), mode = "undirected", diag = FALSE,weighted=TRUE)

color<-adjustcolor(vertexcolor[-which(rowSums(adjacency)==0)], alpha.f = 1)
adjustcolor("grey", alpha.f = 0.3)
color[which(vertexcolor[-which(rowSums(adjacency)==0)]%in%c("purple"))]<-adjustcolor(vertexcolor[-which(rowSums(adjacency)==0)][which(vertexcolor[-which(rowSums(adjacency)==0)]%in%c("purple"))], alpha.f = 0.8)

color[which(vertexcolor[-which(rowSums(adjacency)==0)]%in%c("blue"))]<-adjustcolor(vertexcolor[-which(rowSums(adjacency)==0)][which(vertexcolor[-which(rowSums(adjacency)==0)]%in%c("blue"))], alpha.f = 0.6)
pdf("/home/gridsan/jyun/network/overlap_network_10_15_p_4_1_3.pdf", wi = 10, he = 14)
set.seed(5)
plot(STABLEgraph1,vertex.size= c(size/150+2)[-which(rowSums(adjacency)==0)],vertex.label.cex=0.7,vertex.label.family="Helvetica",vertex.label.color="black",vertex.color= color, edge.color ="white", vertex.frame.color =adjustcolor(vertexcolor[-which(rowSums(adjacency)==0)], alpha.f = 0.9),mark.groups = group_ids,mark.border.width=0.2,mark.border = group_color,mark.col = group_color_FULL,mark.expand=2)

legend("bottomright", title="node (module type)    polygon (overlapping type)",ncol = 2, c(" Bd21c", " Bd21d"," Bd31c"," Bd31d","","conserved","genotypic","environmental","GxE (3)","GxE (1)"),
    pch = c(21,21,21,21,21, 22,22,22,22,22),
    col = c("white","white","white","white","white","white","white","white","white","white"), pt.bg = c("#0000FF99","#FF0000FF","#FFA500FF","#A020F0CC","white","#26C6DA","#FFFF00","#4CAF50","#EC407A","#FFAB91"), pt.cex = 1,border=NA,cex=0.8)

dev.off();



###Module individual GO, 
###group common GO, genes, and annotations
#Module traits, DGE

#install.packages("rbioapi",repos="https://cran.r-project.org")
library(rbioapi)
library(backports)

GO1<-matrix("0",length(levels(as.factor(bd21cmembership[,4]))[-1]),200)
k=1
for (i in levels(as.factor(bd21cmembership[,4]))[-1]){
a1<-a2<-a3<-a4<-a1_terms<-a2_terms <-a3_terms <-a4_terms <-NULL
gene<-colnames(matrixbd)[which(bd21cmembership[,4]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a1<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a1_terms<-a1[which(a1$plus_minus=="+"),]$term.id
}
gene<-colnames(matrixbd)[which(bd21cmembership[,2]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a2<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a2_terms<-a2[which(a2$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a2_terms, a1_terms))
}
gene<-colnames(matrixbd)[which(bd21cmembership[,7]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a3<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a3_terms<-a3[which(a3$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a3_terms, a1_terms))
}

if(i %in% bd21cmembership_mid_breakdown_list_first){
list<-as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])
list<-list[-which(list==i)]
	for (j in 1:length(list)){
		gene<-colnames(matrixbd)[which(bd21cmembership[,2]==list[j])]
		a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a4<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a4_terms<-a4[which(a4$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a4_terms, a1_terms))
}
		
}
}
if (length(a1_terms)>0){
	GO1[k,1:length(a1_terms)]<-a1_terms

}
	k=k+1
}



GO2<-matrix("0",length(levels(as.factor(bd21dmembership[,4]))[-1]),200)
k=1
for (i in levels(as.factor(bd21dmembership[,4]))[-1]){
a1<-a2<-a3<-a4<-a1_terms<-a2_terms <-a3_terms <-a4_terms <-NULL
gene<-colnames(matrixbd)[which(bd21dmembership[,4]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a1<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a1_terms<-a1[which(a1$plus_minus=="+"),]$term.id
}
gene<-colnames(matrixbd)[which(bd21dmembership[,2]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a2<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a2_terms<-a2[which(a2$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a2_terms, a1_terms))
}
gene<-colnames(matrixbd)[which(bd21dmembership[,7]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a3<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a3_terms<-a3[which(a3$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a3_terms, a1_terms))
}

if(i %in% bd21dmembership_mid_breakdown_list_first){
list<-as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])
list<-list[-which(list==i)]
	for (j in 1:length(list)){
		gene<-colnames(matrixbd)[which(bd21dmembership[,2]==list[j])]
		a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a4<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a4_terms<-a4[which(a4$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a4_terms, a1_terms))
}
		
}
}
if (length(a1_terms)>0){
	GO2[k,1:length(a1_terms)]<-a1_terms

}
	k=k+1
}



GO3<-matrix("0",length(levels(as.factor(bd31cmembership[,4]))[-1]),200)
k=1
for (i in levels(as.factor(bd31cmembership[,4]))[-1]){
a1<-a2<-a3<-a4<-a1_terms<-a2_terms <-a3_terms <-a4_terms <-NULL
gene<-colnames(matrixbd)[which(bd31cmembership[,4]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a1<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a1_terms<-a1[which(a1$plus_minus=="+"),]$term.id
}
gene<-colnames(matrixbd)[which(bd31cmembership[,2]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a2<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a2_terms<-a2[which(a2$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a2_terms, a1_terms))
}
gene<-colnames(matrixbd)[which(bd31cmembership[,7]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a3<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a3_terms<-a3[which(a3$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a3_terms, a1_terms))
}

if(i %in% bd31cmembership_mid_breakdown_list_first){
list<-as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list ==i)], split = ",")[[1]])
list<-list[-which(list==i)]
	for (j in 1:length(list)){
		gene<-colnames(matrixbd)[which(bd31cmembership[,2]==list[j])]
		a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a4<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a4_terms<-a4[which(a4$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a4_terms, a1_terms))
}
		
}
}
if (length(a1_terms)>0){
	GO3[k,1:length(a1_terms)]<-a1_terms

}
	k=k+1
}
GO4<-matrix("0",length(levels(as.factor(bd31dmembership[,4]))[-1]),200)
k=1
for (i in levels(as.factor(bd31dmembership[,4]))[-1]){
a1<-a2<-a3<-a4<-a1_terms<-a2_terms <-a3_terms <-a4_terms <-NULL
gene<-colnames(matrixbd)[which(bd31dmembership[,4]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a1<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a1_terms<-a1[which(a1$plus_minus=="+"),]$term.id
}
gene<-colnames(matrixbd)[which(bd31dmembership[,2]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a2<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a2_terms<-a2[which(a2$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a2_terms, a1_terms))
}
gene<-colnames(matrixbd)[which(bd31dmembership[,7]==i)]

a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a3<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a3_terms<-a3[which(a3$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a3_terms, a1_terms))
}

if(i %in% bd31dmembership_mid_breakdown_list_first){
list<-as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])
list<-list[-which(list==i)]
	for (j in 1:length(list)){
		gene<-colnames(matrixbd)[which(bd31dmembership[,2]==list[j])]
		a<-substr(gene,6,nchar(gene)-2)
BRADI<-rep("BRADI_", length(a))
v3<-rep("v3", length(a))
gene <- paste0(BRADI,a,v3)
if (length(gene)>0){
	a4<-rba_panther_enrich(genes=gene,organism=15368,annot_dataset='GO:0003674',test_type='FISHER',correction='FDR',cutoff = 0.1)$result
a4_terms<-a4[which(a4$plus_minus=="+"),]$term.id
a1_terms<-unique(c(a4_terms, a1_terms))
}
		
}
}
if (length(a1_terms)>0){
	GO4[k,1:length(a1_terms)]<-a1_terms

}
	k=k+1
}


library(GO.db)
GO<-rbind(GO1,GO2,GO3,GO4)



write.csv(GO, file=paste("/home/gridsan/jyun/network/GO.csv",sep=""))
GO <-read.csv("/home/gridsan/jyun/network/GO.csv", header=T,row.names=1)


GOlower<-golabel<-matrix("0",dim(GO)[1],200)
rownames(GO)<-rownames(GOlower)<-rownames(golabel)<-colnames(adjacency)
for (i in 1:dim(GO)[1]){
	if(length(which(GO[i,]!="0"))>0){
	k=1
for (j in 1: length(which(GO[i,]!="0"))){
	if(GO[i,j]%in%names(as.list(GOBPOFFSPRING)) ){
if ("TRUE"%in%(GO[i,]%in%c(as.character(GOBPOFFSPRING [GO[i,j]])))){
	}else{
	GOlower[i,k]<-GO[i,j]
	k=k+1
}
}
}
}
}
for (i in 1:dim(GO)[1]){
for (j in 1:length(which(GOlower[i,]!="0"))){
	if(length(which(names(Term(GOTERM))==GOlower[i,j]))>0){
golabel[i,j]<-Term(GOTERM)[which(names(Term(GOTERM))==GOlower[i,j])]
}
}
}
golabel


write.csv(golabel, file=paste("/home/gridsan/jyun/network/golabel.csv",sep=""))
golabel <-read.csv("/home/gridsan/jyun/network/golabel.csv", header=T,row.names=1)

chloroplastgenes<-colnames(matrixbd)[which(substr (colnames(matrixbd),1,12)=="lcl.LT558597")]
chloroplastgenes

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,4] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
	if(i%in%bd21cmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
	if(i %in%bd21cmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
j=1
if(i%in% bd21cmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,7] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast3<-ifelse(p<0.0005,1,0)

chloroplast_bd21c<-ifelse(chloroplast1+ chloroplast2+ chloroplast2_1+ chloroplast2_2+ chloroplast2_3+chloroplast3>1,1,0)


chloroplastgenes<-colnames(matrixbd)[which(substr (colnames(matrixbd),1,12)=="lcl.LT558597")]
chloroplastgenes

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,4] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
	if(i%in%bd21dmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
	if(i %in%bd21dmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
j=1
if(i%in% bd21dmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,7] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast3<-ifelse(p<0.0005,1,0)

chloroplast_bd21d<-ifelse(chloroplast1+ chloroplast2+ chloroplast2_1+ chloroplast2_2+ chloroplast2_3+chloroplast3>1,1,0)



chloroplastgenes<-colnames(matrixbd)[which(substr (colnames(matrixbd),1,12)=="lcl.LT558597")]
chloroplastgenes

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,4] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
	if(i%in%bd31cmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
	if(i %in%bd31cmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
j=1
if(i%in% bd31cmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,7] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast3<-ifelse(p<0.0005,1,0)

chloroplast_bd31c<-ifelse(chloroplast1+ chloroplast2+ chloroplast2_1+ chloroplast2_2+ chloroplast2_3+chloroplast3>1,1,0)

chloroplastgenes<-colnames(matrixbd)[which(substr (colnames(matrixbd),1,12)=="lcl.LT558597")]
chloroplastgenes

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,4] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
	if(i%in%bd31dmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
	if(i %in%bd31dmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
j=1
if(i%in% bd31dmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,7] ==i)]%in%chloroplastgenes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%chloroplastgenes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
chloroplast3<-ifelse(p<0.0005,1,0)

chloroplast_bd31d<-ifelse(chloroplast1+ chloroplast2+ chloroplast2_1+ chloroplast2_2+ chloroplast2_3+chloroplast3>1,1,0)



chloroplast<-c(chloroplast_bd21c, chloroplast_bd21d, chloroplast_bd31c, chloroplast_bd31d)

golabel[,dim(golabel)[2]]= ifelse(chloroplast==1,"chloroplast_genes",0)

rownames(GO)<-colnames(adjacency)


adjacency_prepare <-colnames(finalmatirx_1_number_keep)
commongeneannotation_sort<-commongeneannotation_unsort <-matrix(0,length(adjacency_prepare),700)
commongene_percentage<-matrix(0,length(adjacency_prepare),4)
commongene<-matrix(0,length(adjacency_prepare),700)


for (i in 1:length(adjacency_prepare)){
	a=b=c=d="10"
	if(strsplit (adjacency_prepare[i], split = " ")[[1]][1]!="-"){
		a=colnames(matrixbd)[which(bd21cmembership[,4]==strsplit (adjacency_prepare[i], split = " ")[[1]][1])]
		if(length(a)==0){
		a=colnames(matrixbd)[which(bd21cmembership[,2]==strsplit (adjacency_prepare[i], split = " ")[[1]][1])]			
		}
			}
		if(strsplit (adjacency_prepare[i], split = " ")[[1]][2]!="-"){
		b=colnames(matrixbd)[which(bd21dmembership[,4]==strsplit (adjacency_prepare[i], split = " ")[[1]][2])]
		if(length(b)==0){
		b=colnames(matrixbd)[which(bd21dmembership[,2]==strsplit (adjacency_prepare[i], split = " ")[[1]][2])]			
		}
			}
		if(strsplit (adjacency_prepare[i], split = " ")[[1]][3]!="-"){
		c=colnames(matrixbd)[which(bd31cmembership[,4]==strsplit (adjacency_prepare[i], split = " ")[[1]][3])]
			if(length(c)==0){
		c=colnames(matrixbd)[which(bd31cmembership[,2]==strsplit (adjacency_prepare[i], split = " ")[[1]][3])]			
		}
		}
		if(strsplit (adjacency_prepare[i], split = " ")[[1]][4]!="-"){
		d=colnames(matrixbd)[which(bd31dmembership[,4]==strsplit (adjacency_prepare[i], split = " ")[[1]][4])]
			if(length(d)==0){
		d=colnames(matrixbd)[which(bd31dmembership[,2]==strsplit (adjacency_prepare[i], split = " ")[[1]][4])]			
		}
		}
	if("10"%in%c(a,b,c,d)){
		commongene[i,1:length(plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,2]==c(4-plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,1]=="10"),2])),1])]<-as.vector(plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,2]==c(4-plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,1]=="10"),2])),1])
		l<-plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,2]==c(4-plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,1]=="10"),2])),1]
		}else{
commongene[i,1:length(plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,2]==4),1])]<-as.vector(plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,2]==4),1])
			l<-plyr::count(c(a,b,c,d))[which(plyr::count(c(a,b,c,d))[,2]==4),1]
	}
	

	h<-rep(1,length(l))
if(length(l)>0){
RNA_Brachy_ANNOTATION3<-read.csv("/home/gridsan/jyun/network/annotation_info_PART3.csv", header=T)
#~/Desktop/codes jwafs/annotation_info_PART3.csv
for (j in 1: length(l)){
if(substr(as.character(l[j]),1,12)%in% RNA_Brachy_ANNOTATION3[,1]){
h[j]<-as.character(RNA_Brachy_ANNOTATION3[which(RNA_Brachy_ANNOTATION3[,1]==substr(l[j],1,12))[1],3])

h[j]<-ifelse("lcl.LT558597"==substr(l[j],1,12),"chloroplast gene",h[j])
}
}
}
commongeneannotation_sort[i,1:length(l)]<-sort(h)
commongeneannotation_unsort[i,1:length(l)]<-h

	
commongene_percentage[i,]<-c(ifelse(length(a)>1,length(l)/length(a),0),ifelse(length(b)>1,length(l)/length(b),0),ifelse(length(c)>1,length(l)/length(c),0),ifelse(length(d)>1,length(l)/length(d),0))
}




GOcommonlower <-GOcommonlabel<-GOcommon<-matrix(0,dim(finalmatirx_1_number_keep)[2],100)
for (i in 1:length(colnames(finalmatirx_1_number_keep))){
	names<-NULL
	if(strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[1]!="-"){
		names<-c(names,paste0("Bd21c-",strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][1])	)
	}
		if(strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[2]!="-"){
		names<-c(names,paste0("Bd21d-",strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][2])	)
	}
		if(strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[3]!="-"){
		names<-c(names,paste0("Bd31c-",strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][3])	)
	}
		if(strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[4]!="-"){
		names<-c(names,paste0("Bd31d-",strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]][4])	)
	}
	length<-length(plyr::count(as.vector(as.matrix(GO[which(rownames(GO)%in%names),])))[which(plyr::count(as.vector(as.matrix(GO[which(rownames(GO)%in%names),])))[,2]==sum(ifelse(strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]]=="-",0,1))),1])
	if(length>0){
	GOcommon[i,1: length]<-as.character(plyr::count(as.vector(as.matrix(GO[which(rownames(GO)%in%names),])))[which(plyr::count(as.vector(as.matrix(GO[which(rownames(GO)%in%names),])))[,2]==sum(ifelse(strsplit (colnames(finalmatirx_1_number_keep)[i], split = " ")[[1]]=="-",0,1))),1])
}
}
for (i in 1:dim(GOcommon)[1]){
	if(length(which(GOcommon[i,]!="0"))>0){
	k=1
for (j in which(GOcommon[i,]!="0")){
	if(GOcommon[i,j]%in%names(as.list(GOBPOFFSPRING)) ){
if ("TRUE"%in%(GOcommon[i,]%in%c(as.character(GOBPOFFSPRING [GOcommon[i,j]])))){
	}else{
	GOcommonlower[i,k]<-GOcommon[i,j]
	k=k+1
}
}
}
}
}



for (i in 1:dim(GOcommonlower)[1]){
for (j in which(GOcommonlower[i,]!="0")){
	if(length(which(names(Term(GOTERM))== GOcommonlower[i,j]))>0){
GOcommonlabel[i,j]<-Term(GOTERM)[which(names(Term(GOTERM))== GOcommonlower[i,j])]
print(i)

}
if(GOcommonlower[i,j]=="chloroplast_genes"){
	GOcommonlabel[i,j]<-GOcommonlower[i,j]
}
}
}
GOcommonlabel

write.csv(cbind(colnames(finalmatirx_1_number_keep),GOcommonlabel), file=paste("/home/gridsan/jyun/network/overlapping_GOcommonlabel.csv",sep=""))
write.csv(cbind(colnames(finalmatirx_1_number_keep),commongene_percentage), file=paste("/home/gridsan/jyun/network/overlapping_commongene_percentage.csv",sep=""))
write.csv(cbind(colnames(finalmatirx_1_number_keep),commongene), file=paste("/home/gridsan/jyun/network/overlapping_commongene.csv",sep=""))
write.csv(cbind(colnames(finalmatirx_1_number_keep),commongeneannotation_sort), file=paste("/home/gridsan/jyun/network/overlapping_commongeneannotation_sort.csv",sep=""))
write.csv(cbind(colnames(finalmatirx_1_number_keep),commongeneannotation_unsort), file=paste("/home/gridsan/jyun/network/overlapping_commongeneannotation_unsort.csv",sep=""))
write.csv(cbind(oldnames,golabel), file=paste("/home/gridsan/jyun/network/go_label.csv",sep=""))
Data_otherinfo_updatedbatch<-read.csv("/home/gridsan/jyun/network/Data_otherinfo_updatedbatch.csv", header=T,row.names=1)
trait_data<-Data_otherinfo_updatedbatch[match(rownames(matrixbd),paste0("X",Data_otherinfo_updatedbatch$Experiment.location.ID)),c(15:18,20,22:24,26:30)]



trait_p_all_1<-matrix(0,length(levels(as.factor(bd21cmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_all_1<-matrix(0,length(levels(as.factor(bd21cmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[,which(bd21cmembership[,4]==i)])$x[,1]	

	for (j in 1:dim(trait_data)[2]){
	trait_p_all_1[k,j]<-c(-log(cor.test(eigengene, trait_data[,j])$p.value,10))
	trait_cor_all_1[k,j]<-cor.test(eigengene, trait_data[,j])$estimate 
}
	k=k+1
}
trait_p_all_2<-matrix(0,length(levels(as.factor(bd21dmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_all_2<-matrix(0,length(levels(as.factor(bd21dmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[,which(bd21dmembership[,4]==i)])$x[,1]	
	for (j in 1:dim(trait_data)[2]){
	trait_p_all_2[k,j]<-c(-log(cor.test(eigengene, trait_data[,j])$p.value,10))
	trait_cor_all_2[k,j]<-cor.test(eigengene, trait_data[,j])$estimate 
}
	k=k+1
}

trait_p_all_3<-matrix(0,length(levels(as.factor(bd31cmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_all_3<-matrix(0,length(levels(as.factor(bd31cmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[,which(bd31cmembership[,4]==i)])$x[,1]	
	for (j in 1:dim(trait_data)[2]){
	trait_p_all_3[k,j]<-c(-log(cor.test(eigengene, trait_data[,j])$p.value,10))
	trait_cor_all_3[k,j]<-cor.test(eigengene, trait_data[,j])$estimate 
}
	k=k+1
}

trait_p_all_4<-matrix(0,length(levels(as.factor(bd31dmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_all_4<-matrix(0,length(levels(as.factor(bd31dmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[,which(bd31dmembership[,4]==i)])$x[,1]	
	for (j in 1:dim(trait_data)[2]){
	trait_p_all_4[k,j]<-c(-log(cor.test(eigengene, trait_data[,j])$p.value,10))
	trait_cor_all_4[k,j]<-cor.test(eigengene, trait_data[,j])$estimate 
}
	k=k+1
}
trait_p_all<-rbind(trait_p_all_1,trait_p_all_2,trait_p_all_3,trait_p_all_4)
trait_cor_all<-rbind(trait_cor_all_1,trait_cor_all_2,trait_cor_all_3,trait_cor_all_4)


trait_p_librarywise_1<-matrix(0,length(levels(as.factor(bd21cmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_librarywise_1<-matrix(0,length(levels(as.factor(bd21cmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[1:48,which(bd21cmembership[,4]==i)])$x[,1]	
	for (j in 1:dim(trait_data)[2]){
	trait_p_librarywise_1[k,j]<-c(-log(cor.test(eigengene, trait_data[1:48,j])$p.value,10))
	trait_cor_librarywise_1[k,j]<-cor.test(eigengene, trait_data[1:48,j])$estimate 
}
	k=k+1
}
trait_p_librarywise_2<-matrix(0,length(levels(as.factor(bd21dmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_librarywise_2<-matrix(0,length(levels(as.factor(bd21dmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[49:96,which(bd21dmembership[,4]==i)])$x[,1]	
	for (j in 1:dim(trait_data)[2]){
	trait_p_librarywise_2[k,j]<-c(-log(cor.test(eigengene, trait_data[49:96,j])$p.value,10))
	trait_cor_librarywise_2[k,j]<-cor.test(eigengene, trait_data[49:96,j])$estimate 
}
	k=k+1
}

trait_p_librarywise_3<-matrix(0,length(levels(as.factor(bd31cmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_librarywise_3<-matrix(0,length(levels(as.factor(bd31cmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[97:144,which(bd31cmembership[,4]==i)])$x[,1]	
	for (j in 1:dim(trait_data)[2]){
	trait_p_librarywise_3[k,j]<-c(-log(cor.test(eigengene, trait_data[97:144,j])$p.value,10))
	trait_cor_librarywise_3[k,j]<-cor.test(eigengene, trait_data[97:144,j])$estimate 
}
	k=k+1
}

trait_p_librarywise_4<-matrix(0,length(levels(as.factor(bd31dmembership[,4]))[-1]),dim(trait_data)[2])
trait_cor_librarywise_4<-matrix(0,length(levels(as.factor(bd31dmembership[,4]))[-1]),dim(trait_data)[2])
k=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
	eigengene <-prcomp(matrixbd[145:192,which(bd31dmembership[,4]==i)])$x[,1]	
	for (j in 1:dim(trait_data)[2]){
	trait_p_librarywise_4[k,j]<-c(-log(cor.test(eigengene, trait_data[145:192,j])$p.value,10))
	trait_cor_librarywise_4[k,j]<-cor.test(eigengene, trait_data[145:192,j])$estimate 
}
	k=k+1
}
trait_p_librarywise<-rbind(trait_p_librarywise_1,trait_p_librarywise_2,trait_p_librarywise_3,trait_p_librarywise_4)
trait_cor_librarywise<-rbind(trait_cor_librarywise_1,trait_cor_librarywise_2,trait_cor_librarywise_3,trait_cor_librarywise_4)

textMatrix_library= matrix(paste(signif(trait_cor_librarywise, 1), "\n(",
                        signif(trait_p_librarywise, 1), ")", sep = ""),nrow=dim(trait_cor_librarywise)[1])
textMatrix_all = matrix(paste(signif(trait_cor_all, 1), "\n(",
                        signif(trait_p_all, 1), ")", sep = ""),nrow=dim(trait_cor_librarywise)[1])
                        
colnames(trait_p_librarywise)<-colnames(trait_cor_librarywise)<-colnames(textMatrix_library)<-colnames(textMatrix_all)<-colnames(trait_p_all)<-colnames(trait_cor_all)<-c("day5 water usage","Stomatal conductance (gsw)","Net photosynthesis rate", "Transpiration rate", "Leaf water content", "Starch", "High D.P. Fructan","Protein", "Glucose","Fructose","Sucrose", "Low D.P. Fructan","Total amino acids")
rownames(trait_p_librarywise)<-rownames(trait_cor_librarywise)<-rownames(textMatrix_library)<-rownames(textMatrix_all)<-rownames(trait_p_all)<-rownames(trait_cor_all)<-c(paste0("Bd21c_",levels(as.factor(bd21cmembership[,4]))[-1]),paste0("Bd21d_",levels(as.factor(bd21dmembership[,4]))[-1]),paste0("Bd3-1c_",levels(as.factor(bd31cmembership[,4]))[-1]),paste0("Bd3-1d_",levels(as.factor(bd31dmembership[,4]))[-1]))
                        
                        

library(gplots)
my_palette <- colorRampPalette(c("red", "yellow"))(n = 299)

pdf("/home/gridsan/jyun/network/Module_trait_relationships_in_all_samples_12_15.pdf", wi = 11, he = 10)
par(mar = c(14, 8.8, 3, 2.2));
heatmap.2(trait_p_all,
  cellnote = textMatrix_all,  # same data set for cell labels
  main = "Module--trait relationships in all samples", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none", 
  notecex = 0.2,        # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
 Rowv = FALSE,
    Colv = FALSE, 
    srtCol=30)            # turn off column clustering
dev.off();

pdf("/home/gridsan/jyun/network/Module_trait_relationships_in_library_samples12_15.pdf", wi = 11, he = 10)
par(mar = c(14, 8.8, 3, 2.2));
heatmap.2(trait_p_librarywise,
  cellnote = textMatrix_library,  # same data set for cell labels
  main = "Module--trait relationships in specific library", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none", 
  notecex = 0.2,        # turns off trace lines inside the heat map
  margins =c(12,9),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier    # enable color transition at specified limits
  dendrogram="none",     # only draw a row dendrogram
 Rowv = FALSE,
    Colv = FALSE, 
    srtCol=30)            # turn off column clustering
dev.off();


write.csv(textMatrix_library, file="/home/gridsan/jyun/network/trait_module_Matrix_library_correlation_-logp.csv")

##for DGE in any parameter settings
#g
diffg1111 <-read.csv("/home/gridsan/jyun/network/g_extended.csv", header=T)
e1111genes<-diffe1111[which(diffe1111[,7]<0.01),1]
#gxe
diff_ge <-read.csv("/home/gridsan/jyun/network/gxe_extended.csv", header=T)
e1111genes <-diff_ge[,1]
#e
diffe1111 <-read.csv("/home/gridsan/jyun/network/e_extended.csv", header=T)
e1111genes<-diffe1111[which(diffe1111[,7]<0.01),1]


c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,4] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11111<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
	if(i%in%bd21cmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
	if(i %in%bd21cmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
j=1
if(i%in% bd21cmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==as.numeric(strsplit(bd21cmembership_mid_breakdown_list[which(bd21cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,7] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11113<-ifelse(p<0.0005,1,0)

e1111_bd21c<-ifelse(e11111+ e11112+ e11112_1+ e11112_2+ e11112_3+e11113>1,1,0)


c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,4] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11111<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
	if(i%in%bd21dmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
	if(i %in%bd21dmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
j=1
if(i%in% bd21dmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==as.numeric(strsplit(bd21dmembership_mid_breakdown_list[which(bd21dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,7] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11113<-ifelse(p<0.0005,1,0)

e1111_bd21d<-ifelse(e11111+ e11112+ e11112_1+ e11112_2+ e11112_3+e11113>1,1,0)



c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,4] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11111<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
	if(i%in%bd31cmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
	if(i %in%bd31cmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
j=1
if(i%in% bd31cmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==as.numeric(strsplit(bd31cmembership_mid_breakdown_list[which(bd31cmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,7] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11113<-ifelse(p<0.0005,1,0)

e1111_bd31c<-ifelse(e11111+ e11112+ e11112_1+ e11112_2+ e11112_3+e11113>1,1,0)


c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,4] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,4] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11111<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
	if(i%in%bd31dmembership_mid_breakdown_list_first){
		
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[1])])

		
	}
	
		j=j+1	
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_1<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
	if(i %in%bd31dmembership_mid_breakdown_list_first){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[2])])
}
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_2<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
j=1
if(i%in% bd31dmembership_mid_breakdown_list_first){
	if(length(as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]]))==3){
		c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==as.numeric(strsplit(bd31dmembership_mid_breakdown_list[which(bd31dmembership_mid_breakdown_list_first ==i)], split = ",")[[1]])[3])])
}
j=j+1
}
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112_3<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,7] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,7] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11113<-ifelse(p<0.0005,1,0)

e1111_bd31d<-ifelse(e11111+ e11112+ e11112_1+ e11112_2+ e11112_3+e11113>1,1,0)

e1111<-c(e1111_bd21c, e1111_bd21d, e1111_bd31c, e1111_bd31d)

##for DGE in most conserved settings:
#g
diffg1111 <-read.csv("/home/gridsan/jyun/network/g_extended.csv", header=T)
e1111genes<-diffe1111[which(diffe1111[,7]<0.01),1]
#gxe
diff_ge <-read.csv("/home/gridsan/jyun/network/gxe_extended.csv", header=T)
e1111genes <-diff_ge[,1]
#e
diffe1111 <-read.csv("/home/gridsan/jyun/network/e_extended.csv", header=T)
e1111genes<-diffe1111[which(diffe1111[,7]<0.01),1]

c<-matrix(0,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21cmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21cmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)
e1111_bd21c<-e11112

c<-matrix(0,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd21dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd21dmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd21dmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd21dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd21dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd21dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)


e1111_bd21d<-e11112

c<-matrix(0,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31cmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31cmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31cmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31cmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31cmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31cmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)
e1111_bd31c<-e11112

c<-matrix(0,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,2)
j=1
for (i in as.numeric(levels(as.factor(bd31dmembership[,4]))[-1])){
c[j,1]<-length(which(colnames(matrixbd)[which(bd31dmembership[,2] ==i)]%in%e1111genes))
c[j,2]<-length(colnames(matrixbd)[which(bd31dmembership[,2] ==i)])
j=j+1
}
c[(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1,]<-c(length(which(colnames(matrixbd)%in%e1111genes)),16207)
p<-rep(0,(length(levels(as.factor(bd31dmembership[,4]))[-1])))
for (i in 1:(length(levels(as.factor(bd31dmembership[,4]))[-1]))){
p[i]<-fisher.test(c[c(i,(length(levels(as.factor(bd31dmembership[,4]))[-1]))+1),], alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value
}
e11112<-ifelse(p<0.0005,1,0)

e1111_bd31d<-e11112

c(e1111_bd21c, e1111_bd21d, e1111_bd31c, e1111_bd31d)


