##individual module data
##DAG generation
##edge DGE changes, edge, 
##plots and summarization


##get data of the module
	moduellist<-NULL
for(i in levels(as.factor(bd21cmembership[,2]))[-1]){
	if((length(which(bd21cmembership[,2]==i))<70)&(length(which((colnames(matrixbd)%in%ge)&(bd21cmembership[,2]==i)))>0)){
		print(i)
		print(length(which((colnames(matrixbd)%in%ge)&(bd21cmembership[,2]==i))))
	print(length(which(bd21cmembership[,2]==i)))
	moduellist<-c(moduellist,paste0("bd21c","_",i))
	}
	}
RNA_Brachy_other <-read.csv("/home/gridsan/jyun/dci/Data_otherinfo_updatedbatch.csv", header=T)
#RNA_Brachy_other <-read.csv("~/Desktop/codes\ jwafs/Data_otherinfo_updatedbatch.csv", header=T)
RNA_Brachy_other<-RNA_Brachy_other[-210,c(2,3,4,16,17,18,19,21,23,24,25,27,28,29,30)]

RNA_Brachy_other$Accessiontreatment<-paste0(RNA_Brachy_other$Accession, RNA_Brachy_other$Treatment)
RNA_Brachy_other$Accessiontreatment<-as.factor(RNA_Brachy_other$Accessiontreatment)


datExpr <-read.csv("~/Desktop/all things/expected_count_with_cloroplast_genes_normalized.csv", header=T,row.names=1)

for (i in 1:length(moduellist)){
	if(strsplit(moduellist[i],"_")[[1]][1]=="bd21c"){
		a=strsplit(moduellist[i],"_")[[1]][2]
write.csv(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD21normal"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21cmembership[,2]==a)])],paste0("/home/gridsan/jyun/dci/","bd21c_",a,"_bd21c.csv"))
write.csv(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD21drought"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21cmembership[,2]==a)])],paste0("/home/gridsan/jyun/dci/","bd21c_",a,"_bd21d.csv"))
write.csv(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD3-1normal"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21cmembership[,2]==a)])],paste0("/home/gridsan/jyun/dci/","bd21c_",a,"_bd31c.csv"))
write.csv(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD3-1drought"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21cmembership[,2]==a)])],paste0("/home/gridsan/jyun/dci/","bd21c_",a,"_bd31d.csv"))	

	}
	}

##add glucose Bd31d8, Bd31d12, Bd21d13, Bd21d12
glucose<-log(RNA_Brachy_other$glucose.nmol.Glc.equivalents.mg)
water<-RNA_Brachy_other$day5_waterusage

moduellist="bd21d_13"
		if(strsplit(moduellist,"_")[[1]][1]=="bd21d"){
		a=strsplit(moduellist,"_")[[1]][2]
		b<-na.omit(cbind(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD21normal"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21dmembership[,2]==a)])], glucose[which(RNA_Brachy_other$Accessiontreatment=="BD21normal")]))
			colnames(b)[dim(b)[2]]<-"glucose"
write.csv(b,paste0("/home/gridsan/jyun/dci/","bd21d_",a,"_glucose_real_bd21c.csv"))
b<-na.omit(cbind(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD21drought"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21dmembership[,2]==a)])], glucose[which(RNA_Brachy_other$Accessiontreatment=="BD21drought")]))
	colnames(b)[dim(b)[2]]<-"glucose"
write.csv(b,paste0("/home/gridsan/jyun/dci/","bd21d_",a,"_glucose_real_bd21d.csv"))
b<-na.omit(cbind(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD3-1normal"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21dmembership[,2]==a)])], glucose[which(RNA_Brachy_other$Accessiontreatment=="BD3-1normal")]))
	colnames(b)[dim(b)[2]]<-"glucose"
write.csv(b,paste0("/home/gridsan/jyun/dci/","bd21d_",a,"_glucose_real_bd31c.csv"))
b<-na.omit(cbind(datExpr[which(RNA_Brachy_other$Accessiontreatment=="BD3-1drought"),which(colnames(datExpr)%in%colnames(matrixbd)[which(bd21dmembership[,2]==a)])], glucose[which(RNA_Brachy_other$Accessiontreatment=="BD3-1drought")]))
	colnames(b)[dim(b)[2]]<-"glucose"
write.csv(b,paste0("/home/gridsan/jyun/dci/","bd21d_",a,"_glucose_real_bd31d.csv"))		}


##the DAG is run from CausalAnalysis-Bd31d12_example.ipynb

##
moduellist<-c( "bd21c_6" , "bd21c_7" , "bd21c_8" , "bd21c_9" , "bd21c_11" ,"bd21c_13",
"bd21c_14", "bd21c_17", "bd21c_21", "bd21c_22" ,"bd21c_23" ,"bd21d_7", 
 "bd21d_12", "bd21d_13" ,"bd21d_16","bd31d_7_real" , "bd31d_8_real" ,
"bd31d_9_real" , "bd31d_11_real" ,"bd31d_12_real" ,"bd31d_14_real", "bd31d_18_real",
 "bd31d_19_real"  ,"bd31c_8"  ,"bd31c_13" ,"bd31c_14","bd21c_10","bd31d_14_trait","bd21d_2","bd21d_12_glucose_real","bd21d_13_glucose_real","bd31d_8_glucose_real","bd31d_12_glucose_real","bd31c_10","bd31c_16","bd31d_18_water_real","bd31d_11_water_real")
matrixbd<-read.csv("~/Desktop/draft_SI/draft_data/matrixbd_for_groupwiseWGCNA.csv", header=T,row.names=1)
RNA_Brachy_other <-read.csv("~/Desktop/codes\ jwafs/Data_otherinfo_updatedbatch.csv", header=T,row.names=1)

glucose<-log(RNA_Brachy_other$glucose.nmol.Glc.equivalents.mg)
water<-RNA_Brachy_other$day5_waterusage

matrixbd<-cbind(matrixbd[,1:16207],glucose[match(rownames(matrixbd),paste0("X",RNA_Brachy_other[,1]))])
matrixbd<-cbind(matrixbd,water[match(rownames(matrixbd),paste0("X",RNA_Brachy_other[,1]))])
 colnames(matrixbd)[16208]<-"glucose"
 colnames(matrixbd)[16209]<-"water"
 
 n=which(moduellist=="bd31c_16")

 
 DGEall<-matrix(0,6,37)
bigtable<-matrix(0,130,130)

for(n in 1:length(moduellist)){

DGECHANGESUMMARY<-NULL
NNN=1
 table<-NULL

	#input data
	filename<-moduellist[n]
	genes2 <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"_bd31d.csv"), header=T, row.names = 1)

	pdf(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"TEST10_withnotraits_anno_multi_loose_onedi_test3.pdf"), width=12, height = 12)
moduellist[n]
genes2 <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"_bd21d.csv"), header=T, row.names = 1)
genes <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"_bd21c.csv"), header=T, row.names = 1)
gegenes <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/adjacency_matrix_", filename,".csv"), header=F)
rownames(gegenes)<-colnames(gegenes)<-colnames(genes)
gegenes<-t(gegenes)
pquantile=ifelse(dim(gegenes)[2]>40,0.97,0.95)
cutoff=0.5
adjMat=as.matrix(ifelse(gegenes>=min(cutoff,round(quantile(as.matrix(gegenes), pquantile),3)),1,0))

rownames(adjMat)<-colnames(adjMat)<-colnames(genes)
genes1 <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"_bd21c.csv"), header=T, row.names = 1)
genes1<-genes1[,which(colnames(genes1)%in%colnames(genes))]
genes2 <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"_bd21d.csv"), header=T, row.names = 1)
genes2<-genes2[,which(colnames(genes2)%in%colnames(genes))]
genes3 <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"_bd31c.csv"), header=T, row.names = 1)
genes3<-genes3[,which(colnames(genes3)%in%colnames(genes))]
genes4 <-read.csv(paste0("~/Desktop/chapter1_code_data_updated_12_21/dag/", filename,"_bd31d.csv"), header=T, row.names = 1)
genes4<-genes4[,which(colnames(genes4)%in%colnames(genes))]

#two complete same genes in one module
if(n==2){
	genes1<-genes1[,-54]
		genes2<-genes2[,-54]
			genes3<-genes3[,-54]
				genes4<-genes4[,-54]
				gegenes<-gegenes[-54,-54]
}
#annotation
genes<-rbind(genes1, genes2, genes3, genes4)
normgenes<-as.data.frame(rbind(scale(genes1, center=TRUE, scale=TRUE),scale(genes2, center=TRUE, scale= TRUE),scale(genes3, center=TRUE, scale= TRUE),scale(genes4, center=TRUE, scale= TRUE)))
l<-colnames(genes)
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
cbind(l,h)
#DGE

groupe<-c(rep(1,dim(genes1)[1]),rep(0,dim(genes2)[1]),rep(1,dim(genes3)[1]),rep(0,dim(genes4)[1]))
groupg<-c(rep(1,dim(genes1)[1]+dim(genes2)[1]),rep(0,dim(genes3)[1]+dim(genes4)[1]))
diff_ge <-read.csv("~/Desktop/draft_SI/draft_data/gxe_extended.csv", head=TRUE)
ge <-diff_ge[,1]
diff_g <-read.csv("~/Desktop/draft_SI/draft_data/g_extended.csv", head=TRUE)
g<-diff_g[,1]
diff_g <-read.csv("~/Desktop/draft_SI/draft_data/e_extended.csv", head=TRUE)
e<-diff_g[,1]

genes_g_e<-colnames(genes)[which((colnames(genes)%in%g)& (colnames(genes)%in%e)&(!(colnames(genes)%in%ge)))]
genes_g<-colnames(genes)[which((colnames(genes)%in%g)&(!(colnames(genes)%in%e))&(!(colnames(genes)%in%ge)))]
genes_e<-colnames(genes)[which((colnames(genes)%in%e)& (!(colnames(genes)%in%g))&(!(colnames(genes)%in%ge)))]
genes_ge<-colnames(genes)[which(colnames(genes)%in%ge)]
#network
library(igraph)
library(Rgraphviz)
am.graph<-new("graphAM", adjMat= as.matrix(adjMat), edgemode="directed")
#color nodes according to DGE
eAttrs <- list()
nAttrs<- list()
Cololnames<-eAttrs$fillcolor <- NULL
#nAttrs$fillcolor=c(rep("goldenrod2",length(which(colnames(genes)%in% genes_g_e))),rep("aquamarine2",length(which(colnames(genes)%in% genes_g))),rep("deepskyblue2",length(which(colnames(genes)%in% genes_e))),rep("coral2",length(which(colnames(genes)%in% genes_ge))))
#names(nAttrs$fillcolor)<-c(colnames(adjMat)[which(colnames(genes)%in% genes_g_e)],colnames(adjMat)[which(colnames(genes)%in% genes_g)],colnames(adjMat)[which(colnames(genes)%in% genes_e)],colnames(adjMat)[which(colnames(genes)%in% genes_ge)])
nAttrs$fillcolor=c(rep("yellow",length(which(colnames(genes)%in% genes_g_e))),rep("green",length(which(colnames(genes)%in% genes_g))),rep("blue",length(which(colnames(genes)%in% genes_e))),rep("red",length(which(colnames(genes)%in% genes_ge))))
names(nAttrs$fillcolor)<-c(colnames(adjMat)[which(colnames(genes)%in% genes_g_e)],colnames(adjMat)[which(colnames(genes)%in% genes_g)],colnames(adjMat)[which(colnames(genes)%in% genes_e)],colnames(adjMat)[which(colnames(genes)%in% genes_ge)])
names_eAttrs_fillcolor <-edge_dge<-rep("black",dim(which(adjMat>0,arr.ind=TRUE))[1])
##ase

for(i in 1:length(substr(colnames(genes1),1,12))){
ase<-rep("-",dim(genes)[2])
for(i in 1:dim(genes)[2]){
	if(substr(colnames(genes),1,12)[i]%in%substr(rownames(data_ase),1,12)){
	ase[i]<-"NA"
	}

	if(substr(colnames(genes),1,12)[i]%in%substr(rownames(data_ase),1,12)){
	ase[i]<-paste0(data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),1], data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),2], data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),3], data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),4], data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),5], data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),6], data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),7], data_ase[which(substr(rownames(data_ase),1,12)==substr(colnames(genes),1,12)[i]),8])
	}
	
}
#plot ASE
names1<-Names<-substr(colnames(genes1),1,12)[i]
if(length(which(rownames(normalized_counts_orignal)== Names))>0){
data<-cbind(sample_real_id_matrix2$Accession, sample_real_id_matrix2 $Treatment, sample_real_id_matrix2 $generation, c(normalized_counts_orignal[which(rownames(normalized_counts_orignal)== Names),]))
colnames(data)<-c("Accession","Treatment","generation","expression")
write.csv(data,"~/Data_temperatl.csv")
data <-read.csv("~/Data_temperatl.csv", header=T)
boxplot(expression ~ Treatment*Accession*generation,data=data,xlab=paste0(Names ,"-",cis[which(rownames(normalized_counts_orignal)== Names)],cisxe[which(rownames(normalized_counts_orignal)== Names)],cis_both[which(rownames(normalized_counts_orignal)== Names)],cis_e[which(rownames(normalized_counts_orignal)== Names)],trans[which(rownames(normalized_counts_orignal)== Names)],transxe[which(rownames(normalized_counts_orignal)== Names)],trans_both[which(rownames(normalized_counts_orignal)== Names)],trans_e[which(rownames(normalized_counts_orignal)== Names)]))

if(length(which(substr(colnames(matrixbd),1,12)== Names))>0){

plot(matrixbd[,which(substr(colnames(matrixbd),1,12)== Names)])
}
}
}

#check edges 
diag(gegenes)<-0
a<-b<-NULL

ppp=0.05
r2=0.5
gegenes0<-gegenes
	groupg0<-groupg
	groupe0<-groupe
	genes0<-genes
#if need subsample, ended up not using it, if use can also change pp to be more times
SIZE=45

interceptchangegene<-rep("",dim(gegenes)[1])
	changeedge<-changeedge_label<-lis <-NULL
for (i in 1:dim(gegenes)[1]){
	#
		#print(colnames(genes1)[i])
	if(length(which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile))))>0){
for (j in 1:length(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))])){
print(paste0(colnames(genes1)[i],"~",colnames(genes1)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))
#print(car::Anova(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes),0.94)))[j]]), collapse="+"),sep = ""),data= genes), type=2))

blist <-matrix(1,2,3)
for(pp in 1:2){

			groupg<-groupg0
	groupe<-groupe0
	genes <-genes0
	subsample<-c(sample(which((groupg0==1)&(groupe0==1)),SIZE,replace=TRUE),sample(which((groupg0==1)&(groupe0==0)), SIZE,replace= TRUE),sample(which((groupg0==0)&(groupe0==1)), SIZE,replace= TRUE),sample(which((groupg0==0)&(groupe0==0)), SIZE,replace= TRUE))

		#	groupg<-groupg[subsample]
#groupe<-groupe[subsample]
#	genes <-genes[subsample,]
#cook distance to remove outlier samples
	# mod<-lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes)
# cooksd_rem <- which(cooks.distance(mod)>4/192)
# if(length(cooksd_rem)>0){
		# groupg<-groupg[-cooksd_rem]
	# groupe<-groupe[-cooksd_rem]
	# genes <-genes[-cooksd_rem,]
	# }else{
		# groupg<-groupg
	# groupe<-groupe
	# genes <-genes			
		# }
	result<-car::Anova(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes), type=2)	
#blist is the slope variation p-values
	b<-as.vector(result[c(which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupg")),paste0(c("groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe")),paste0(c("groupe:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe:groupg")),paste0(c("groupe:groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j])))),4])
blist[pp,]<-b
blist[is.na(blist)]<-1
}

#interceptchangegene is the intercept variation type with p=0.05

	interceptchangegene[i]<-paste0(ifelse((result[which(rownames(result)=="groupe:groupg"),4]<0.05)&(b[3]>0.05),"gxe","-"),ifelse((result[which(rownames(result)=="groupe"),4]>0.05)&(result[which(rownames(result)=="groupg"),4]<0.05)&(result[which(rownames(result)=="groupe:groupg"),4]>0.05)&(b[3]>0.05)&(b[1]>0.05),"g","-"),ifelse((result[which(rownames(result)=="groupe"),4]<0.05)&(result[which(rownames(result)=="groupg"),4]>0.05)&(result[which(rownames(result)=="groupe:groupg"),4]>0.05)&(b[3]>0.05)&(b[2]>0.05),"e","-"))
##rsquare
	op<-signif(summary(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes))$r.squared,digits=2)

##for each pairs, look more closely the regression and plot 
##check first multiple regression 

	m=which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))[j]
	k=i	
	
					groupg<-groupg0
	groupe<-groupe0
	genes <-genes0
				result<-car::Anova(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes), type=2)

resulti<-summary(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes))$coefficients

slopep <-signif(-log(as.vector(result[c(which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupg")),paste0(c("groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe")),paste0(c("groupe:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe:groupg")),paste0(c("groupe:groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j])))),4]),10),digits=2)

slopei <-signif(as.vector(resulti[c(which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupg")),paste0(c("groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe")),paste0(c("groupe:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe:groupg")),paste0(c("groupe:groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j])))),1]),digits=1)
	
intercep <-signif(-log(as.vector(result[c(which(rownames(as.matrix(result))%in%c("groupg","groupe","groupe:groupg"))),4])[c(2,1,3)],10),digits = 2)
intercei <-signif(as.vector(resulti[c(which(rownames(as.matrix(resulti))%in%c("groupg","groupe","groupe:groupg"))),1])[c(2,1,3)],digits = 1)
condition <-c(rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48),rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48))
value <-c(genes[,k],genes[,m])
Gene1<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[m])]
Gene2<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[k])]

genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))

value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value)
Data<-as.data.frame(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
Data<-na.omit(Data)
dgecp<-signif(-log(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4],10),digits = 2)
dgeci<-signif(summary(lm(value~condition*gene*genotype,data= Data))$coefficients[c(5,7,8),1],digits = 1)	


	##slope,r2, intercept,dge,dge change, ase	
		mgene<-paste0(colnames(adjMat)[m],"-int", intercei[1],"p", intercep[1],"-", intercei[2],"p", intercep[2],"-", intercei[3],"p", intercep[3] ,ifelse(ase[m]=="-","-","ase"),ase[m])
	kgene<-paste0(colnames(adjMat)[k],"-m-sl","-", slopei[1],"p", slopep[1],"-", slopei[2],"p", slopep[2],"-", slopei[3],"p", slopep[3] , ifelse(ase[k]=="-","-","ase"),ase[k])
		listdata0<-c(mgene,kgene,slopep[1],slopep[2],slopep[3])
		
			plot(genes[which((groupe==1)&(groupg==1)),m],genes[which((groupe==1)&(groupg==1)),k],col="blue",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene,main=paste0(paste0(ifelse((slopep[1]>1.29)&(op>0.5),"g","-"),"~",ifelse((slopep[2]>1.29)&(op>0.5),"e","-"),"~",ifelse((slopep[3]>1.29)&(op>0.5),"gxe","-"),"~","dge", dgeci[1],"p",dgecp[1],"-",dgeci[2],"p",dgecp[2],"-",dgeci[3],"p",dgecp[3] ),"r",op))
	#print(type)
	abline(lm(genes[which((groupe==1)&(groupg==1)),k]~genes[which((groupe==1)&(groupg==1)),m]),col="blue")
	par(new=TRUE)
plot(genes[which((groupe==0)&(groupg==1)),m],genes[which((groupe==0)&(groupg==1)),k],col="red",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==0)&(groupg==1)),k]~genes[which((groupe==0)&(groupg==1)),m]),col="red")
	par(new=TRUE)
	plot(genes[which((groupe==1)&(groupg==0)),m],genes[which((groupe==1)&(groupg==0)),k],col="purple",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==1)&(groupg==0)),k]~genes[which((groupe==1)&(groupg==0)),m]),col="purple")
	par(new=TRUE)
	plot(genes[which((groupe==0)&(groupg==0)),m],genes[which((groupe==0)&(groupg==0)),k],col="orange",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==0)&(groupg==0)),k]~genes[which((groupe==0)&(groupg==0)),m]),col="orange")
	
##single regression 
m=which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))[j]
	k=i	
					groupg<-groupg0
	groupe<-groupe0
	genes <-genes0
		op<-signif(summary(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]), collapse="+"),sep = ""),data= genes))$r.squared,digits=2)
	
result<-car::Anova(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]), collapse="+"),sep = ""),data= genes), type=2)

resulti<-summary(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]), collapse="+"),sep = ""),data= genes))$coefficients
slopep <-signif(-log(as.vector(result[c(which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupg")),paste0(c("groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe")),paste0(c("groupe:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(result))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe:groupg")),paste0(c("groupe:groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j])))),4]),10),digits=2)




slopei <-signif(as.vector(resulti[c(which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupg")),paste0(c("groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe")),paste0(c("groupe:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe:groupg")),paste0(c("groupe:groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j])))),1]),digits=1)


	
intercep <-signif(-log(as.vector(result[c(which(rownames(as.matrix(result))%in%c("groupg","groupe","groupe:groupg"))),4])[c(2,1,3)],10),digits = 2)
intercei <-signif(as.vector(resulti[c(which(rownames(as.matrix(resulti))%in%c("groupg","groupe","groupe:groupg"))),1])[c(2,1,3)],digits = 1)
condition <-c(rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48),rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48))
value <-c(genes[,k],genes[,m])
Gene1<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[m])]
Gene2<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[k])]

genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))

value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value)

Data<-as.data.frame(Data)
Data<-na.omit(Data)

write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
Data<-na.omit(Data)
dgecp<-signif(-log(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4],10),digits = 2)
dgeci<-signif(summary(lm(value~condition*gene*genotype,data= Data))$coefficients[c(5,7,8),1],digits = 1)	

	##slope,r2, intercept,dge,dge change, ase	
		mgene<-paste0(colnames(adjMat)[m],"-int", intercei[1],"p", intercep[1],"-", intercei[2],"p", intercep[2],"-", intercei[3],"p", intercep[3] ,"-",ase[m])
	kgene<-paste0(colnames(adjMat)[k],"-m-sl","-", slopei[1],"p", slopep[1],"-", slopei[2],"p", slopep[2],"-", slopei[3],"p", slopep[3] , ase[k],"-")
			plot(genes[which((groupe==1)&(groupg==1)),m],genes[which((groupe==1)&(groupg==1)),k],col="blue",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene,main=paste0(paste0("dge", dgeci[1],"p",dgecp[1],"-",dgeci[2],"p",dgecp[2],"-",dgeci[3],"p",dgecp[3] ),"r",op))
	#print(type)
	abline(lm(genes[which((groupe==1)&(groupg==1)),k]~genes[which((groupe==1)&(groupg==1)),m]),col="blue")
	par(new=TRUE)
plot(genes[which((groupe==0)&(groupg==1)),m],genes[which((groupe==0)&(groupg==1)),k],col="red",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==0)&(groupg==1)),k]~genes[which((groupe==0)&(groupg==1)),m]),col="red")
	par(new=TRUE)
	plot(genes[which((groupe==1)&(groupg==0)),m],genes[which((groupe==1)&(groupg==0)),k],col="purple",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==1)&(groupg==0)),k]~genes[which((groupe==1)&(groupg==0)),m]),col="purple")
	par(new=TRUE)
	plot(genes[which((groupe==0)&(groupg==0)),m],genes[which((groupe==0)&(groupg==0)),k],col="orange",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==0)&(groupg==0)),k]~genes[which((groupe==0)&(groupg==0)),m]),col="orange")
	#print(c(mgene,kgene,slopep[1],slopep[2],slopep[3]))
	listdata1<-c(mgene,kgene,slopep[1],slopep[2],slopep[3])

###correlation by scale gene expression first

m=which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))[j]
	k=i	
					groupg<-groupg0
	groupe<-groupe0
	genes <-genes0
		op<-signif(summary(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes))$r.squared,digits=2)
	gene1a<-scale(genes[,k])
		gene1b<-scale(genes[,m])
		genes1all<-cbind(gene1a, gene1b, groupe, groupg)
		colnames(genes1all)[1:2]<-c(colnames(genes)[k],colnames(genes)[m])
		result<-car::Anova(lm(paste0(colnames(genes)[k],"~groupe*groupg*",colnames(genes)[m]),data= as.data.frame(genes1all)), type=2)

genesscaled<-as.data.frame(rbind(scale(genes[1:dim(genes1)[1],]),scale(genes[(dim(genes1)[1]+1):(dim(genes1)[1]+dim(genes2)[1]),]),scale(genes[(dim(genes1)[1]+dim(genes2)[1]+1):(dim(genes1)[1]+dim(genes2)[1]+dim(genes3)[1]),]),scale(genes[(dim(genes1)[1]+dim(genes2)[1]+dim(genes3)[1]+1):(dim(genes1)[1]+dim(genes2)[1]+dim(genes3)[1]+dim(genes4)[1]),])))
resulti<-summary(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]), collapse="+"),sep = ""),data= genesscaled))$coefficients
slopep <-signif(-log(as.vector(resulti[c(which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupg")),paste0(c("groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe")),paste0(c("groupe:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe:groupg")),paste0(c("groupe:groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j])))),4]),10),digits=2)





slopei <-signif(as.vector(resulti[c(which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupg")),paste0(c("groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe")),paste0(c("groupe:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))),which(rownames(as.matrix(resulti))%in%c(paste0(colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j],c(":groupe:groupg")),paste0(c("groupe:groupg:"),colnames(genes)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j])))),1]),digits=1)


	
intercep <-signif(-log(as.vector(result[c(which(rownames(as.matrix(result))%in%c("groupg","groupe","groupe:groupg"))),4])[c(2,1,3)],10),digits = 2)
intercei <-signif(as.vector(resulti[c(which(rownames(as.matrix(resulti))%in%c("groupg","groupe","groupe:groupg"))),1])[c(2,1,3)],digits = 1)
condition <-c(rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48),rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48))
value <-c(genes[,k],genes[,m])
Gene1<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[m])]
Gene2<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[k])]

genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))

value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value)
Data<-as.data.frame(Data)
Data<-na.omit(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
Data<-na.omit(Data)
dgecp<-signif(-log(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4],10),digits = 2)
dgeci<-signif(summary(lm(value~condition*gene*genotype,data= Data))$coefficients[c(5,7,8),1],digits = 1)	


	##slope,r2, intercept,dge,dge change, ase	
		mgene<-paste0(colnames(adjMat)[m],"-int", intercei[1],"p", intercep[1],"-", intercei[2],"p", intercep[2],"-", intercei[3],"p", intercep[3] ,"-",ase[m])
	kgene<-paste0(colnames(adjMat)[k],"-cor","-", slopei[1],"p", slopep[1],"-", slopei[2],"p", slopep[2],"-", slopei[3],"p", slopep[3] , ase[k],"-")
			plot(genes[which((groupe==1)&(groupg==1)),m],genes[which((groupe==1)&(groupg==1)),k],col="blue",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene,main=paste0(paste0("dge", dgeci[1],"p",dgecp[1],"-",dgeci[2],"p",dgecp[2],"-",dgeci[3],"p",dgecp[3] ),"r",op))
	#print(type)
	
	abline(lm(genes[which((groupe==1)&(groupg==1)),k]~genes[which((groupe==1)&(groupg==1)),m]),col="blue")
	par(new=TRUE)
plot(genes[which((groupe==0)&(groupg==1)),m],genes[which((groupe==0)&(groupg==1)),k],col="red",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==0)&(groupg==1)),k]~genes[which((groupe==0)&(groupg==1)),m]),col="red")
	par(new=TRUE)
	plot(genes[which((groupe==1)&(groupg==0)),m],genes[which((groupe==1)&(groupg==0)),k],col="purple",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==1)&(groupg==0)),k]~genes[which((groupe==1)&(groupg==0)),m]),col="purple")
	par(new=TRUE)
	plot(genes[which((groupe==0)&(groupg==0)),m],genes[which((groupe==0)&(groupg==0)),k],col="orange",ylim=c(min(genes[,k]),max(genes[,k])),xlim=c(min(genes[,m]),max(genes[,m])),xlab= mgene,ylab= kgene)
	abline(lm(genes[which((groupe==0)&(groupg==0)),k]~genes[which((groupe==0)&(groupg==0)),m]),col="orange")
	
###differential DGE	
	condition <-c(rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48),rep("Drought",48), rep("Control",48), rep("Drought",48), rep("Control",48))
value <-c(genes[,k],genes[,m])
Gene1<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[m])]
Gene2<-matrixbd[,which(colnames(matrixbd)==colnames(adjMat)[k])]

genotype <-c(rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48),rep("Bd21",48), rep("Bd21",48), rep("Bd3-1",48), rep("Bd3-1",48))
gene <-c(rep("Gene1",192), rep("Gene2", 192))

value<-c(Gene1,Gene2)
Data <-cbind(condition,gene, genotype, value)
Data<-as.data.frame(Data)
Data<-na.omit(Data)
write.csv(Data,"~/Data_temperatl.csv")
Data<-read.csv("~/Data_temperatl.csv", header=T,row.names=1)
Data<-na.omit(Data)
dgecp<-car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4]
dgeci<-signif(summary(lm(value~condition*gene*genotype,data= Data))$coefficients[c(5,7,8),1],digits = 1)	
print(c(dgecp, dgeci))

if(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4][3]<0.001){
	
	DGECHANGESUMMARY[NNN]<-"GxE"
}else{
	if((car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4][2]<0.001)&(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4][1]<0.001)){
			DGECHANGESUMMARY[NNN]<-"G+E"
	}else{
		if(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4][1]<0.001){
	
	DGECHANGESUMMARY[NNN]<-"E"
}else{
	if(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4][2]<0.001){
	
	DGECHANGESUMMARY[NNN]<-"G"
}
if(car::Anova(lm(value~condition*gene*genotype,data= Data), type=2)[c(4,6,7),4][2]>=0.001){
	
	DGECHANGESUMMARY[NNN]<-"non"
}
	}
	}
}
NNN=NNN+1

DGEall[1,n]<-filename
DGEall[2,n]<-ifelse("G"%in%plyr::count(DGECHANGESUMMARY)[,1],plyr::count(DGECHANGESUMMARY)[which(plyr::count(DGECHANGESUMMARY)[,1]=="G"),2],0)
DGEall[3,n]<-ifelse("E"%in%plyr::count(DGECHANGESUMMARY)[,1],plyr::count(DGECHANGESUMMARY)[which(plyr::count(DGECHANGESUMMARY)[,1]=="E"),2],0)
DGEall[4,n]<-ifelse("G+E"%in%plyr::count(DGECHANGESUMMARY)[,1],plyr::count(DGECHANGESUMMARY)[which(plyr::count(DGECHANGESUMMARY)[,1]=="G+E"),2],0)
DGEall[5,n]<-ifelse("GxE"%in%plyr::count(DGECHANGESUMMARY)[,1],plyr::count(DGECHANGESUMMARY)[which(plyr::count(DGECHANGESUMMARY)[,1]=="GxE"),2],0)
DGEall[6,n]<-length(DGECHANGESUMMARY)



#GET BACK TO THE SUMMARY DATA
groupg<-groupg0
	groupe<-groupe0
	genes <-genes0
		op<-signif(summary(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes))$r.squared,digits=2)
		result<-car::Anova(lm(paste(colnames(genes)[i],"~",paste(paste0("groupe*groupg*",colnames(genes)[which(gegenes[,i]>min(cutoff,quantile(as.matrix(gegenes), pquantile)))]), collapse="+"),sep = ""),data= genes), type=2)
		
lis <-rbind(lis,colSums(ifelse(blist <ppp,1,0)))	#if (min(b)<p){
	changeedge<-c(changeedge, paste0(colnames(genes1)[i],"~",colnames(genes1)[which(gegenes[,i]> min(cutoff,quantile(as.matrix(gegenes), pquantile)))][j]))
	#changeedge_label<-c(changeedge_label,paste0(ifelse(colSums(ifelse(blist <ppp,1,0))[3]>60,"gxe",""), ifelse((colSums(ifelse(blist <ppp,1,0))[1]>60)&(colSums(ifelse(blist <ppp,1,0))[2]>60)&(colSums(ifelse(blist <ppp,1,0))[3]<60),"g+e",""),ifelse((colSums(ifelse(blist <ppp,1,0))[1]>60)&(colSums(ifelse(blist <ppp,1,0))[2]<60)&(colSums(ifelse(blist <ppp,1,0))[3]<60),"g",""),ifelse((colSums(ifelse(blist <ppp,1,0))[1]<60)&(colSums(ifelse(blist <ppp,1,0))[2]>60)&(colSums(ifelse(blist <ppp,1,0))[3]<60),"e","")))
		
		changeedge_label<-c(changeedge_label,paste0(ifelse(((colSums(ifelse(blist <ppp,1,0))[3]>1)&(op> r2)),"gxe",""), ifelse(((colSums(ifelse(blist <ppp,1,0))[1]>1)&(colSums(ifelse(blist <ppp,1,0))[2]>1)&(colSums(ifelse(blist <ppp,1,0))[3]<1)&(op> r2)),"g+e",""),ifelse(((colSums(ifelse(blist <ppp,1,0))[1]>1)&(colSums(ifelse(blist <ppp,1,0))[2]<1)&(colSums(ifelse(blist <ppp,1,0))[3]<1)&(op> r2)),"g",""),ifelse(((colSums(ifelse(blist <ppp,1,0))[1]<1)&(colSums(ifelse(blist <ppp,1,0))[2]>1)&(colSums(ifelse(blist <ppp,1,0))[3]<1)&(op> r2)),"e","")))
}
	}	
}

##COLLECT DATA FOR BIGTABLE
table<-c(table,length(which(changeedge_label=="e"))/length(changeedge_label))
table<-c(table,length(which(changeedge_label=="g"))/length(changeedge_label))
table<-c(table,length(which(changeedge_label=="g+e"))/length(changeedge_label))
table<-c(table,length(which(changeedge_label=="gxe"))/length(changeedge_label))

if((length(which(colnames(genes)%in% genes_ge))>2)&(length(which(colnames(genes)%in% genes_ge))<length(colnames(genes))-2)){
	#print(t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% ge), which(colnames(genes)%in% ge)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% ge), which(colnames(genes)%in% ge)])>0,1,0),alternative ="greater"))

table<-c(table,t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% genes_ge), which(colnames(genes)%in% genes_ge)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% genes_ge), which(colnames(genes)%in% genes_ge)])>0,1,0),alternative ="greater")$p.value)

}else{
table<-c(table,"na")	
	}
	
table<-c(table,length(which(colnames(genes)%in% genes_ge)),length(which(interceptchangegene[which(colnames(genes)%in% genes_ge)]=="gxe--"))/length(interceptchangegene[which(colnames(genes)%in% genes_ge)]),fisher.test(rbind(c(length(which(interceptchangegene[which(colnames(genes)%in% genes_ge)]=="gxe--")),length(interceptchangegene[which(colnames(genes)%in% genes_ge)])),
c(length(which(interceptchangegene[-which(colnames(genes)%in% genes_ge)]=="gxe--")),length(interceptchangegene[-which(colnames(genes)%in% genes_ge)])))
, alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value

)
plyr::count(interceptchangegene[which(colnames(genes)%in% genes_ge)])
#print(interceptchangegene)
#print(length(which(colnames(genes)%in% ge))/length(genes))

if((length(which(colnames(genes)%in% genes_g))>2)&(length(which(colnames(genes)%in% genes_g))<length(colnames(genes))-2)){
	#print(t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% g), which(colnames(genes)%in% g)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% g), which(colnames(genes)%in% g)])>0,1,0),alternative ="greater"))
	table<-c(table,t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% genes_g), which(colnames(genes)%in% genes_g)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% genes_g), which(colnames(genes)%in% genes_g)])>0,1,0),alternative ="greater")$p.value)

}else{
table<-c(table,"na")	
	}
table<-c(table,length(which(colnames(genes)%in% genes_g)),length(which(interceptchangegene[which(colnames(genes)%in% genes_g)]=="-g-"))/length(interceptchangegene[which(colnames(genes)%in% genes_g)]),fisher.test(rbind(c(length(which(interceptchangegene[which(colnames(genes)%in% genes_g)]=="-g-")),length(interceptchangegene[which(colnames(genes)%in% genes_g)])),
c(length(which(interceptchangegene[-which(colnames(genes)%in% genes_g)]=="-g-")),length(interceptchangegene[-which(colnames(genes)%in% genes_g)])))
, alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value)
plyr::count(interceptchangegene[which(colnames(genes)%in% genes_g)])
#print(interceptchangegene)
#print(length(which(colnames(genes)%in% g))/length(genes))
if((length(which(colnames(genes)%in% genes_g_e))>2)&(length(which(colnames(genes)%in% genes_g_e))<length(colnames(genes))-2)){
	#print(t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% g), which(colnames(genes)%in% g)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% g), which(colnames(genes)%in% g)])>0,1,0),alternative ="greater"))
	table<-c(table,t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% genes_g_e), which(colnames(genes)%in% genes_g_e)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% genes_g_e), which(colnames(genes)%in% genes_g_e)])>0,1,0),alternative ="greater")$p.value)

}else{
table<-c(table,"na")	
	}
table<-c(table,length(which(colnames(genes)%in% genes_g_e)),length(which(interceptchangegene[which(colnames(genes)%in% genes_g_e)]=="-g+e-"))/length(interceptchangegene[which(colnames(genes)%in% genes_g_e)]),fisher.test(rbind(c(length(which(interceptchangegene[which(colnames(genes)%in% genes_g_e)]=="-g+e-")),length(interceptchangegene[which(colnames(genes)%in% genes_g_e)])),
c(length(which(interceptchangegene[-which(colnames(genes)%in% genes_g_e)]=="-g+e-")),length(interceptchangegene[-which(colnames(genes)%in% genes_g_e)])))
, alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value)
plyr::count(interceptchangegene[which(colnames(genes)%in% genes_g_e)])
#print(interceptchangegene)
#print(length(which(colnames(genes)%in% g))/length(genes))


if((length(which(colnames(genes)%in% genes_e))>2)&(length(which(colnames(genes)%in% genes_e))<length(colnames(genes))-2)){
	#print(t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% e), which(colnames(genes)%in% e)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% e), which(colnames(genes)%in% e)])>0,1,0),alternative ="greater"))
	
	table<-c(table,t.test(ifelse(rowSums(adjMat[which(colnames(genes)%in% genes_e), which(colnames(genes)%in% genes_e)])>0,1,0),ifelse(rowSums(adjMat[-which(colnames(genes)%in% genes_e), which(colnames(genes)%in% genes_e)])>0,1,0),alternative ="greater")$p.value)

	print(n)
}else{
table<-c(table,"na")	
	}

table<-c(table,length(which(colnames(genes)%in% genes_e)),length(which(interceptchangegene[which(colnames(genes)%in% genes_e)]=="--e"))/length(interceptchangegene[which(colnames(genes)%in% genes_e)]),fisher.test(rbind(c(length(which(interceptchangegene[which(colnames(genes)%in% genes_e)]=="--e")),length(interceptchangegene[which(colnames(genes)%in% genes_e)])),
c(length(which(interceptchangegene[-which(colnames(genes)%in% genes_e)]=="--e")),length(interceptchangegene[-which(colnames(genes)%in% genes_e)])))
, alternative = "greater",
            conf.int = TRUE, conf.level = 0.95,
            simulate.p.value = TRUE, B = 2000)[1]$p.value)

plyr::count(interceptchangegene[which(colnames(genes)%in% genes_e)])
#print(interceptchangegene)
#print(length(which(colnames(genes)%in% e))/length(genes))

bigtable[n,1:length(table)]<- table
print(table)

print(n)
print(plyr::count(DGECHANGESUMMARY))



nAttrs$shape<-rep("box",length(rownames(adjMat)))
names(nAttrs$shape)<-rownames(adjMat)


changeedgeNEW<-NULL
for (i in 1:length(strsplit(changeedge, split = "~"))){
	
	changeedgeNEW<-c(changeedgeNEW,paste0(colnames(adjMat)[which(colnames(genes)==strsplit(changeedge, split = "~")[[i]][1])],"~",colnames(adjMat)[which(colnames(genes)==strsplit(changeedge, split = "~")[[i]][2])]))
}

changeedge1<-NULL
for (i in 1:length(strsplit(changeedgeNEW, split = "~"))){
	
	changeedge1<-c(changeedge1,paste0(strsplit(changeedgeNEW, split = "~")[[i]][2],"~",strsplit(changeedgeNEW, split = "~")[[i]][1]))
}
names_eAttrs_fillcolorNEW<-NULL
for (i in 1:length(strsplit(names_eAttrs_fillcolor, split = "~"))){
	
	names_eAttrs_fillcolorNEW <-c(names_eAttrs_fillcolorNEW,paste0(colnames(adjMat)[which(colnames(genes)==strsplit(names_eAttrs_fillcolor, split = "~")[[i]][1])],"~",colnames(adjMat)[which(colnames(genes)==strsplit(names_eAttrs_fillcolor, split = "~")[[i]][2])]))
}

eAttrs$label<-changeedge_label
eAttrs$label<-DGECHANGESUMMARY
names(eAttrs$label)<-changeedge1

require(Rgraphviz)
Rgraphviz::plot(am.graph, edgeAttrs= eAttrs,nodeAttrs= nAttrs, recipEdges="distinct" )
orgina<-cbind(DGECHANGESUMMARY,names(eAttrs$label))
print(cbind(sort(names(eAttrs$label)), orgina[match(sort(names(eAttrs$label)),orgina[,2]),]))
sort(nAttrs$fillcolor)		
dev.off()

}

write.csv(DGEall,"~/Downloads/DGEall.csv")

write.csv(bigtable[1:37,1:40],"~/Downloads/bigtable.csv")




groupe<-c(rep("c",48),rep("d",48),rep("c",48),rep("d",48))[which(glucose[match(rownames(matrixbd),paste0("X",RNA_Brachy_other$Experiment.location.ID))]>0)]
groupg<-c(rep("bd21",96),rep("bd31",96))[which(glucose[match(rownames(matrixbd),paste0("X",RNA_Brachy_other$Experiment.location.ID))]>0)]
matrixbd[1:5,1:5]

result<-car::Anova(lm(paste("glucose","~",paste(paste0("groupe*groupg*",c("Bradi4g39520.1","Bradi1g22050.1","Bradi1g64687.1","Bradi2g48450.1","Bradi3g51250.1","Bradi4g38380.1","Bradi2g52310.1","Bradi4g42400.1","Bradi5g18470.4")), collapse="+"),sep = ""),data= matrixbd[which(glucose[match(rownames(matrixbd),paste0("X",RNA_Brachy_other$Experiment.location.ID))]>0),]), type=2)	
