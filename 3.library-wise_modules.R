matrixbd<-read.csv("~/network/matrixbd_for_groupwiseWGCNA.csv", header=T,row.names=1)

library(WGCNA)

dim(matrixbd)
###set 

nSets = 5; 
options(stringsAsFactors = FALSE)

# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Bd21n","Bd31d", "Bd31n","Bd21d","all")
shortLabels = c("Bd21n","Bd31d", "Bd31n","Bd21d","all")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr1 = vector(mode = "list", length = nSets)
multiExpr1[[1]] = list(data = na.omit(as.data.frame(matrixbd[1:48,])));
multiExpr1[[2]] = list(data = na.omit(as.data.frame(matrixbd[49:96,])));
multiExpr1[[3]] = list(data = na.omit(as.data.frame(matrixbd[97:144,])));
multiExpr1[[4]] = list(data = na.omit(as.data.frame(matrixbd[c(145:158,160:192),])));
#multiExpr1[[4]] = list(data = na.omit(as.data.frame(matrixbd[c(145:192),])));

multiExpr1[[5]] = list(data = na.omit(as.data.frame(matrixbd)));
names(multiExpr1[[1]]$data) = colnames(matrixbd);
rownames(multiExpr1[[1]]$data) =rownames(na.omit(matrixbd[1:48,]));
names(multiExpr1[[2]]$data) = colnames(matrixbd);
rownames(multiExpr1[[2]]$data) =rownames(na.omit(matrixbd[49:96,]));
names(multiExpr1[[3]]$data) = colnames(matrixbd);
rownames(multiExpr1[[3]]$data) =rownames(na.omit(matrixbd[97:144,]));
names(multiExpr1[[4]]$data) = colnames(matrixbd);
#rownames(multiExpr1[[4]]$data) =rownames(na.omit(matrixbd[c(145:192),]));

rownames(multiExpr1[[4]]$data) =rownames(na.omit(matrixbd[c(145:158,160:192),]));
names(multiExpr1[[5]]$data) = colnames(matrixbd);
rownames(multiExpr1[[5]]$data) =rownames(na.omit(matrixbd));
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr1)
gsg = goodSamplesGenesMS(multiExpr1, verbose = 3);
gsg$allOK



  nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples; 

nSets=4
#
#this is determined following the wgcna code to compare two networks
softPower = 5;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
  adjacencies[set, , ] = abs(cor(multiExpr1[[set]]$data, use = "p"))^softPower;
# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ]);
# Define the reference percentile
scaleP = 0.95 
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                          probs = scaleP, type = 8);
  # Scale the male TOM
if (set>1) { 
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
} 
TOM1<-TOM[1, ,]
Tree1 = hclust(as.dist(1-TOM1), method = "average");
express<-multiExpr1[1]
 unmergedLabelsmid = cutreeDynamic(dendro = Tree1, distM = 1-TOM1,
              deepSplit = 2, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsmid1 = mergeCloseModules(express, unmergedLabelsmid, cutHeight = 0.25, verbose = 3)   
mergedLabelsmid1_c1 = labels2colors(mergedLabelsmid1 $colors)

unmergedLabelsbreak = cutreeDynamic(dendro = Tree1, distM = 1-TOM1,
              deepSplit = 3, cutHeight = 0.984,
              minClusterSize = 15,
              pamStage = FALSE );
mergedLabelsbreak2 = mergeCloseModules(express, unmergedLabelsbreak, cutHeight = 0.25, verbose = 3)
mergedLabelsbreak2_c1 = labels2colors(mergedLabelsbreak2 $colors)

 unmergedLabelsbigger = cutreeDynamic(dendro = Tree1, distM = 1-TOM1,
              deepSplit = 1, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsbigger3 = mergeCloseModules(express, unmergedLabelsbigger, cutHeight = 0.25, verbose = 3)
mergedLabelsbigger3_c1 = labels2colors(mergedLabelsbigger3 $colors)


TOM2<-TOM[2, ,]
Tree2 = hclust(as.dist(1-TOM2), method = "average");
express2<-multiExpr1[2]
 unmergedLabelsmid = cutreeDynamic(dendro = Tree2, distM = 1-TOM2,
              deepSplit = 2, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsmid1 = mergeCloseModules(express2, unmergedLabelsmid, cutHeight = 0.25, verbose = 3)   
mergedLabelsmid1_c2 = labels2colors(mergedLabelsmid1 $colors)

unmergedLabelsbreak = cutreeDynamic(dendro = Tree2, distM = 1-TOM2,
              deepSplit = 3, cutHeight = 0.984,
              minClusterSize = 15,
              pamStage = FALSE );
mergedLabelsbreak2 = mergeCloseModules(express2, unmergedLabelsbreak, cutHeight = 0.25, verbose = 3)
mergedLabelsbreak2_c2 = labels2colors(mergedLabelsbreak2 $colors)

 unmergedLabelsbigger = cutreeDynamic(dendro = Tree2, distM = 1-TOM2,
              deepSplit = 1, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsbigger3 = mergeCloseModules(express2, unmergedLabelsbigger, cutHeight = 0.25, verbose = 3)
mergedLabelsbigger3_c2 = labels2colors(mergedLabelsbigger3 $colors)




TOM3<-TOM[3, ,]
Tree3 = hclust(as.dist(1-TOM3), method = "average");
express3<-multiExpr1[3]
 unmergedLabelsmid = cutreeDynamic(dendro = Tree3, distM = 1-TOM3,
              deepSplit = 2, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsmid1 = mergeCloseModules(express3, unmergedLabelsmid, cutHeight = 0.25, verbose = 3)   
mergedLabelsmid1_c3 = labels2colors(mergedLabelsmid1 $colors)

unmergedLabelsbreak = cutreeDynamic(dendro = Tree3, distM = 1-TOM3,
              deepSplit = 3, cutHeight = 0.984,
              minClusterSize = 15,
              pamStage = FALSE );
mergedLabelsbreak2 = mergeCloseModules(express3, unmergedLabelsbreak, cutHeight = 0.25, verbose = 3)
mergedLabelsbreak2_c3 = labels2colors(mergedLabelsbreak2 $colors)

 unmergedLabelsbigger = cutreeDynamic(dendro = Tree3, distM = 1-TOM3,
              deepSplit = 1, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsbigger3 = mergeCloseModules(express3, unmergedLabelsbigger, cutHeight = 0.25, verbose = 3)
mergedLabelsbigger3_c3 = labels2colors(mergedLabelsbigger3 $colors)







TOM4<-TOM[4, ,]
Tree4 = hclust(as.dist(1-TOM4), method = "average");
express4<-multiExpr1[4]
 unmergedLabelsmid = cutreeDynamic(dendro = Tree4, distM = 1-TOM4,
              deepSplit = 2, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsmid1 = mergeCloseModules(express4, unmergedLabelsmid, cutHeight = 0.25, verbose = 3)   
mergedLabelsmid1_c4 = labels2colors(mergedLabelsmid1 $colors)

unmergedLabelsbreak = cutreeDynamic(dendro = Tree4, distM = 1-TOM4,
              deepSplit = 3, cutHeight = 0.984,
              minClusterSize = 15,
              pamStage = FALSE );
mergedLabelsbreak2 = mergeCloseModules(express4, unmergedLabelsbreak, cutHeight = 0.25, verbose = 3)
mergedLabelsbreak2_c4 = labels2colors(mergedLabelsbreak2 $colors)

 unmergedLabelsbigger = cutreeDynamic(dendro = Tree4, distM = 1-TOM4,
              deepSplit = 1, cutHeight = 0.99,
              minClusterSize = 30,
              pamStage = FALSE );
mergedLabelsbigger3 = mergeCloseModules(express4, unmergedLabelsbigger, cutHeight = 0.25, verbose = 3)
mergedLabelsbigger3_c4 = labels2colors(mergedLabelsbigger3 $colors)

#write.csv( as.vector(mergedLabelsmid1 $colors),"~/network/TOM4_new_smaller_logglucose.csv")
#write.csv( as.vector(mergedLabelsbreak2 $colors),"~/network/TOM4_new_smaller2_logglucose.csv")
#write.csv( as.vector(mergedLabelsbigger3 $colors),"~/network/TOM4_new_smaller3_logglucose.csv")


pdf("/home/gridsan/jyun/network/dendrogramtom4.pdf", wi = 11, he = 16)
plotDendroAndColors(Tree1, cbind(mergedLabelsmid1_c1 , mergedLabelsbreak2_c1, mergedLabelsbigger3_c1), c("middle size", "smaller size","bigger size"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(Tree2, cbind(mergedLabelsmid1_c2 , mergedLabelsbreak2_c2, mergedLabelsbigger3_c2), c("middle size", "smaller size","bigger size"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(Tree3, cbind(mergedLabelsmid1_c3 , mergedLabelsbreak2_c3, mergedLabelsbigger3_c3), c("middle size", "smaller size","bigger size"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(Tree4, cbind(mergedLabelsmid1_c4 , mergedLabelsbreak2_c4, mergedLabelsbigger3_c4), c("middle size", "smaller size","bigger size"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

