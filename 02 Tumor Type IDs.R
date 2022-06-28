# This script uses DESeq2 to calculate the number of DEGs between each samples of one tumor type and all other
# samples. It then uses the boxplot function to find the number of samples in each group for which
# the normalized expression value of each DEG is considered an outlier. It then looks to see if the
# number of identified outliers per sample is an outlier amongst all tumors in that tumor type. 
# I do this in two rounds, starting with the tumor type with the greatest number of DEGs (FS), then PWT and
# finally PNST. This identified one FS, one PWT sample and no PNST as outliers. Next, DEGs are recalculated with 
# the above two outliers removed. To avoid biasing the data in favor of the tumor type with the greatest number of 
# DEGs, the samenumber of DEGs (the most significant 185) are then used to perform unsupervised heirarchical 
# clustering of the samples using the heatmap.2 function. This reassigned the two outliers to different tumor types. 
# Final DEGs are then calculated using the reassigned tumor sample #groups and saved.

library(stringr)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(gplots)
library(plotly)
library(ggpattern)

setwd("C:\\users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data")
setwd("C:\\users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data")

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

mNorm<-function(x){
  mMean<-mean(x)
  mSD<-sd(x)
  (x-mMean)/mSD
}

########  Load the prepared data ##########################

load("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\all_data.rda")
load("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\all_data.rda")

mDat<-all_data
colnames(mDat)
mDat<-mDat[,1:16] #analyze just the cSTS
colnames(mDat)

metaD<-read.csv("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-original.csv")
metaD<-read.csv("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-original.csv")
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

################### find initial DEGs based on just the histological IDs #######
####### load files if already computed #########################################
FS.samples <- rownames(metaD[metaD$Histology.group == "FS",])
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])
l2FC.cutoff <- 1

FS.genes <- GetRes(mDat[,FS.samples],mDat[,c(PNST.samples,PWT.samples)], returnAll = TRUE)
FS.genes.sig<-subset(FS.genes, abs(log2FoldChange) > 1 & padj < .05)

PNST.genes <- GetRes(mDat[,PNST.samples],mDat[,c(FS.samples,PWT.samples)], returnAll = TRUE)
PNST.genes.sig<-subset(PNST.genes, abs(log2FoldChange) > 1 & padj < .05)

PWT.genes <- GetRes(mDat[,PWT.samples],mDat[,c(FS.samples,PNST.samples)], returnAll = TRUE)
PWT.genes.sig<-subset(PWT.genes, abs(log2FoldChange) > 1 & padj < .05)

dds<-DESeqDataSetFromMatrix(mDat, metaD, design = ~1)
mCPM<-rlog(dds)
mCPM<-assay(mCPM, normalized = TRUE)

save(metaD, mDat, mCPM, FS.genes,FS.genes.sig,PNST.genes,PNST.genes.sig,PWT.genes,PWT.genes.sig, file = "C:\\users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\initialDegs.rda")

load("C:\\users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\initialDegs.rda")
load("C:\\users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\initialDegs.rda")
dim(FS.genes.sig) #761 x 6
dim(PNST.genes.sig) #11 x 6
dim(PWT.genes.sig) #233 x 6

FS.samples <- rownames(metaD[metaD$Histology.group == "FS",])
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])

##### find any FS outliers
outlier.matrix<-matrix(0,ncol = length(FS.samples), nrow = nrow(FS.genes.sig))
for (i in 1:nrow(FS.genes.sig)){
  x <- mCPM[rownames(FS.genes.sig)[i],FS.samples]
  OutVals = boxplot(x, plot = FALSE)$out
  OutVals
  outlier.matrix[i,which(x %in% OutVals)] <- 1
}

outlier.sum <- apply(outlier.matrix,2,sum)
rbind(FS.samples,outlier.sum)
outlier.data<-data.frame(IDs = FS.samples, outlierSum = outlier.sum)
write.csv(outlier.data, file = "C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\FS_outlier_sums.csv")
write.csv(outlier.data, file = "C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\FS_outlier_sums.csv")

#get any outliers among the sum of all outlier genes
OutVals = boxplot(outlier.sum, plot = FALSE)$out
OutVals #180
#So, the fifth sample (S2-8, #5) is an outlier.

#### Calculate PWT DEGs, with the FS outlier S2-8 removed
FS.samples <- rownames(metaD[metaD$Histology.group == "FS",]) #5 (S2-8) is the outlier
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])

PWT.genes<-GetRes(mDat[,PWT.samples],mDat[,c(FS.samples[-5],PNST.samples)], returnAll = TRUE)
PWT.genes.sig<-subset(PWT.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(PWT.genes.sig) #now we have significant 521 DEGs, instead of 233

#####PWT outliers
dim(PWT.genes.sig) #521 x 6
outlier.matrix<-matrix(0,ncol = length(PWT.samples), nrow = nrow(PWT.genes.sig))

for (i in 1:nrow(PWT.genes.sig)){
  x <- mCPM[rownames(PWT.genes.sig)[i],PWT.samples]
  OutVals = boxplot(x, plot = FALSE)$out
  OutVals
  outlier.matrix[i,which(x %in% OutVals)] <- 1
}

outlier.sum <- apply(outlier.matrix,2,sum)
outlier.sum
rbind(PWT.samples,outlier.sum)
OutVals = boxplot(outlier.sum, plot = FALSE)$out
OutVals #75
#So, the second sample (S1-10) is an outlier.
outlier.data<-data.frame(IDs = PWT.samples, outlierSum = outlier.sum)

write.csv(outlier.data, file = "C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\PWT_outlier_sums.csv")
write.csv(outlier.data, file = "C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\PWT_outlier_sums.csv")

### recalculate PNST genes after removing the outliers in FS and PWT ##########
FS.samples <- rownames(metaD[metaD$Histology.group == "FS",]) #5 (S2-8) is the outlier
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",]) #2 (S1-10) is the outlier

PNST.genes<-GetRes(mDat[,PNST.samples],mDat[,c(FS.samples[-5],PWT.samples[-2])], returnAll = TRUE)
PNST.genes.sig<-subset(PNST.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(PNST.genes.sig) #now we have significant 185 genes.

#redo analysis of PNST outliers with the 185 significant PNST DEGs
dim(PNST.genes.sig) #185
outlier.matrix<-matrix(0,ncol = length(PNST.samples), nrow = nrow(PNST.genes.sig))
for (i in 1:nrow(PNST.genes.sig)){
  x <- mCPM[rownames(PNST.genes.sig)[i],PNST.samples]
  OutVals = boxplot(x, plot = FALSE)$out
  OutVals
  outlier.matrix[i,which(x %in% OutVals)] <- 1
}
outlier.sum <- apply(outlier.matrix,2,sum)
outlier.sum
outlier.data<-data.frame(IDs = PNST.samples, outlierSum = outlier.sum)
outlier.data
write.csv(outlier.data, file = "C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\PNST_outlier_sums.csv")

OutVals = boxplot(outlier.sum, plot = FALSE)$out
OutVals #NO outliers!!

######################### find DEGs for all tumors with outliers removed ########################

## now, identify DEGs using samples without the two outliers ##
FS.samples <- rownames(metaD[metaD$Histology.group == "FS",])
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])

FS.genes<-GetRes(mDat[,FS.samples[-5]],mDat[,c(PNST.samples,PWT.samples[-2])], returnAll = TRUE)
PWT.genes<-GetRes(mDat[,PWT.samples[-2]],mDat[,c(PNST.samples,FS.samples[-5])], returnAll = TRUE)
PNST.genes<-GetRes(mDat[,PNST.samples],mDat[,c(FS.samples[-5],PWT.samples[-2])], returnAll = TRUE)

FS.genes.sig<-subset(FS.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(FS.genes.sig) #now we have 1714 significant genes.

PWT.genes.sig<-subset(PWT.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(PWT.genes.sig) #now we have 776 significant genes.

PNST.genes.sig<-subset(PNST.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(PNST.genes.sig) #now we have 185 significant genes.

#use heatmap.2 to cluster ALL SAMPLES using the top 185 most significant degs for each sample type.
FS.genes.sig<-FS.genes.sig[order(FS.genes.sig$padj, decreasing = F),]
PWT.genes.sig<-PWT.genes.sig[order(PWT.genes.sig$padj, decreasing = F),]
PNST.genes.sig<-PNST.genes.sig[order(PNST.genes.sig$padj, decreasing = F),]

hGenes<-c(rownames(FS.genes.sig[1:185,]), rownames(PWT.genes.sig[1:185,]), rownames(PNST.genes.sig[1:185,]))
hDat<-mCPM[hGenes,]
colnames(hDat)
names<-colnames(hDat)
names[names %in% FS.samples] <- "FS"
names[names %in% PNST.samples] <- "PNST"
names[names %in% PWT.samples] <- "PWT"
colnames(hDat)<-paste0(colnames(hDat),names)

hDat<-apply(hDat,1,mNorm) #z scores for each gene

colfunc <- colorRampPalette(c("blue", "white", "red"))
myColors<-colfunc(50)

heatmap.2(as.matrix(t(hDat)), Colv = TRUE, Rowv = TRUE, dendrogram = "both", col = myColors, trace = "none", keysize = 1, srtCol=45, key.title = NA, key.ylab = NA, cexRow = 1.25, margins = c(5,15), labRow = NA, colsep = c(4,10), sepwidth = c(.1,.1))

#so, unsupervised heirarchical clustering with the 185 most signficant DEGs for each tumor type
#calculated from data without the above identified outliers 
#assigns S-10(PWT) with the PNST tumors and S-8(FS) with the PWT tumors

###################### redo DESEQ with all tumors (reassigned) to get final tumor-type specific DEGs ############################

metaD<-read.csv("C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-reassigned.csv")
metaD<-read.csv("C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Input Data\\Subject Data\\STS Data for study-ES-reassigned.csv")
dim(metaD)
metaD<-metaD[1:16,1:43]
#View(metaD)
rownames(metaD)<-metaD$Sample.ID
#reorder the metaD data
metaD<-metaD[colnames(mDat),]
all.equal(colnames(mDat),metaD$Sample.ID) #TRUE
table(metaD$Histology.group) #4 FS, 6 PNST, 6 PWT

FS.samples <- rownames(metaD[metaD$Histology.group == "FS",])
PNST.samples <- rownames(metaD[metaD$Histology.group == "PNST",])
PWT.samples <- rownames(metaD[metaD$Histology.group == "PWT",])

PWT.genes <- GetRes(mDat[,PWT.samples],mDat[,c(FS.samples,PNST.samples)], returnAll = TRUE)
PWT.genes.sig<-subset(PWT.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(PWT.genes.sig) #875 as of file produced on 7/20/21 and reported in original submission

FS.genes <- GetRes(mDat[,FS.samples],mDat[,c(PWT.samples,PNST.samples)], returnAll = TRUE)
FS.genes.sig<-subset(FS.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(FS.genes.sig) #2090 as of file produced on 7/20/21, also reported in original submission

PNST.genes <- GetRes(mDat[,PNST.samples],mDat[,c(FS.samples,PWT.samples)], returnAll = TRUE)
PNST.genes.sig<-subset(PNST.genes, abs(log2FoldChange) > 1 & padj < .05)
dim(PNST.genes.sig) #291 as of file produced on 7/20/21, also reported in original submission

writeDir<-"C:\\Users\\Eric Shelden\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric"
writeDir<-"C:\\Users\\idaho\\Dropbox\\Shelden lab\\Canine SARC analysis\\Results\\Tumor type specific gene expression\\Eric"
write.csv(as.data.frame(FS.genes.sig), file = paste0(writeDir,"\\FS_specific_genes_final.csv"))
write.csv(as.data.frame(PNST.genes.sig), file = paste0(writeDir,"\\PNST_specific_genes_final.csv"))
write.csv(as.data.frame(PWT.genes.sig), file = paste0(writeDir,"\\PWT_specific_genes_final.csv"))