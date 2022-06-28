#permutation assay to test goodness of fit for tumor ID assignments. Histology indicated 
#tumor groups of size 5 and 6, and the reassigned tumors had groups of 4 and 6 samples.
#Here, I do 1000 DESEq2 runs using a group of the above size pulled at random from all samples
#then I run DESeq of that group against all other samples. This models, for example, running
#the 5 original FS samples against all other samples in the data set. Then I
#ask how many times the number of resulting DEGs is as great or greater than the number of
#DEGs calculated using the specified groups.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
library(BiocManager)
BiocManager::install("DESeq2")
BiocManager::install("edgeR")

library(stringr)
library(edgeR)
library(DESeq2)
library(ggplot2)

###########################functions##########################################
GetRes<-function(group1, group2, padjCutoff = .05, l2FCCutoff = 1, returnAll = FALSE){
  sIDs<-c(colnames(group1),colnames(group2))
  sGroups<-c(rep("G1",ncol(group1)), rep("G2",ncol(group2)))
  sData<-data.frame(ID = sIDs, group = sGroups)
  dds<-DESeqDataSetFromMatrix(cbind(group1,group2), sData, design = ~ group)
  dds<-DESeq(dds)
  result<-results(dds, contrast = c("group", "G1", "G2"))
  if (returnAll){
    return(result)
  }
  else {
    res_sig<-subset(result, padj < padjCutoff & abs(log2FoldChange) > l2FCCutoff)
    return(res_sig)
  }
}

########  Load the prepared data #########################
baseDir<-"C:\\Users\\idaho\\Dropbox\\"
baseDir<-"C:\\Users\\Eric Shelden\\Dropbox\\"
baseDir<-"C:\\Users\\eshel\\Dropbox\\"
load(str_c(baseDir,"Shelden lab\\Canine SARC analysis\\Input Data\\all_data.rda"))

mFile<-"/home/shelden/Dropbox/Shelden lab/Canine SARC analysis/Input Data/all_data_final.rda"
load(mFile)
mDat<-all_data
colnames(mDat)
mDat<-mDat[,1:16] #just the scrolls
colnames(mDat)


#load the metadata, we are just using the sample names to pull random samples for this analysis.
metaDataFile<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-not-reassigned.csv"
metaDataFile<-"C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-not-reassigned.csv"
metaDataFile<-"C:\\Users\\eshel\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-not-reassigned.csv"
metaDataFile<-"/home/shelden/Dropbox/Shelden lab/Canine SARC analysis/Input Data/Subject Data/STS Data for study-original_edited.csv"
metaD<-read.csv(metaDataFile)
dim(metaD)
metaD<-metaD[1:16,1:43]
#View(metaD)
rownames(metaD)<-metaD$Sample.ID
#reorder the metaD data
metaD<-metaD[colnames(mDat),]
all.equal(colnames(mDat),metaD$Sample.ID) #TRUE
table(metaD$Histology.group)
#unassigned: FS 5, PNST 5, PWT 6
#re-assigned: FS: 4, PNST: 6, PWT: 6

#The number of DEGs found for tumor each tumor type compared
#against the other samples were:
# 761(FS)
# 11(PNST)
# 233(PWT)
# After reassigning the two outliers, The
# number of DEGs was:
# dim(PWT.genes.sig) #875 
# dim(FS.genes.sig) #2090 
# dim(PNST.genes.sig) #291 

#here's the test for all samples
setwd("/home/shelden/Dropbox/Shelden lab/Canine SARC analysis/Results/perm")

nTrials<-1000
pTest<-numeric(nTrials)
#mSamples<-list()
for (i in 1:nTrials){
  sampleSize<-sample(c(5,6),1,replace = F) #get two sample sizes of 5 or 6 at random
  G1<-sample(rownames(metaD),sampleSize, replace = F) #randomly pull out that number of samples
  G2<-setdiff(rownames(metaD),G1)                     #second group is all other samples
  tDegs.sig<-GetRes(mDat[G1],mDat[G2],returnAll = F)
  print(paste(i,dim(tDegs.sig)[1],length(G1)))
  pTest[i]<-dim(tDegs.sig)[1] #112
}
save(pTest, file = "pTest1.rda")

#repeat but use the number of samples in each group after reassignment, ie (4,6)..
nTrials<-1000
pTest<-numeric(nTrials)
#mSamples<-list()
for (i in 1:nTrials){
  sampleSize<-sample(c(4,6),1,replace = F) #get two sample sizes of 5 or 6 at random
  G1<-sample(rownames(metaD),sampleSize, replace = F) #randomly pull out that number of samples
  G2<-setdiff(rownames(metaD),G1)                     #second group is all other samples
  tDegs.sig<-GetRes(mDat[G1],mDat[G2],returnAll = F)
  print(paste(i,dim(tDegs.sig)[1],length(G1)))
  pTest[i]<-dim(tDegs.sig)[1] #112
}
save(pTest, file = "pTest2.rda")

load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\perm\\pTest1.rda")
# dim(FS.genes.sig) #761 x 6
# dim(PNST.genes.sig) #11 x 6
# dim(PWT.genes.sig) #233 x 6
sum(pTest >= 761) #7
sum(pTest >= 11) #332
sum(pTest >= 233) #37
7/1000 #.007
332/1000 #.332
37/1000 #.037
#so as before, the number of significant DEGs is greater than chance for FS and PWT, but not PNST


load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\perm\\pTest2.rda")
# dim(PWT.genes.sig) #875 as of file produced on 7/20/21 and reported in original submission
# dim(FS.genes.sig) #2090 as of file produced on 7/20/21, also reported in original submission
# dim(PNST.genes.sig) #291 as of file produced on 7/20/21, also reported in original submission
sum(pTest >= 2090) #0
sum(pTest >= 875) #2
sum(pTest >= 291) #36
0/10000 #0
2/10000 #.0002
36/1000 #.036
